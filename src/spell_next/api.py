"""spell-next web API and front end.

Run: SPELL_NEXT_DATA=data/build uvicorn spell_next.api:app

Single-gene results are served from the precomputed table; multi-gene searches run
live against the in-memory compendium, a few at a time.
"""
import gzip
import io
import json
import os
import re
import sqlite3
import threading
import warnings
from contextlib import asynccontextmanager
from pathlib import Path

import numpy as np
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, PlainTextResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from spell_next.compendium import Compendium
from spell_next import modconfig
from spell_next.go import ASPECTS, GOIndex
from spell_next.search import search
from spell_next.searchlog import SearchLog

# Resolved once at start-up: a rebuild may repoint ROOT/current, but a running server keeps
# reading the build whose compendium it loaded until it is restarted.
DATA_DIR = Path(os.path.realpath(os.environ.get('SPELL_NEXT_DATA', 'data/build')))
# Monthly search logs; set SPELL_NEXT_LOG_DIR to an empty string to disable.
LOG_DIR = os.environ.get('SPELL_NEXT_LOG_DIR', str(DATA_DIR / 'logs'))
STATIC_DIR = Path(__file__).parent / 'static'
MAX_QUERY_GENES = 50
MAX_LIVE_SEARCHES = 4
NAME_SEP = re.compile(r'[\s,;]+')

state = {}


@asynccontextmanager
async def lifespan(app):
    state['comp'] = Compendium.load(DATA_DIR / 'compendium.npz')
    state['live'] = threading.BoundedSemaphore(MAX_LIVE_SEARCHES)
    state['go'] = load_go_index()
    state['config'] = load_config()
    state['assay_cols'] = assay_columns()
    state['precomputed'] = precomputed()
    state['log'] = SearchLog(LOG_DIR or None)
    yield
    state.clear()


app = FastAPI(title='spell-next', lifespan=lifespan)


def db():
    con = sqlite3.connect('file:%s?mode=ro' % (DATA_DIR / 'spell.sqlite'), uri=True)
    con.row_factory = sqlite3.Row
    return con


def gene_info(con, ids):
    """Return {gene_id: gene dict} for the given SPELL gene IDs."""
    if not ids:
        return {}
    rows = con.execute('SELECT * FROM genes WHERE id IN (%s)' % ','.join('?' * len(ids)), list(ids))
    return {r['id']: dict(r) for r in rows}


def dataset_info(con, ids):
    if not ids:
        return {}
    rows = con.execute('SELECT id, geo_id, pubmed_id, title, first_author, authors, journal, year, tags, num_conds, '
                       'num_genes, species, %s FROM datasets WHERE id IN (%s)'
                       % (state['assay_cols'], ','.join('?' * len(ids))), list(ids))
    return {r['id']: dict(r) for r in rows}


def resolve(con, names):
    """Resolve gene names to SPELL gene IDs: systematic and standard names first, then aliases, then the MOD's
    gene ID (S000001855 or SGD:S000001855)."""
    resolved, unknown, ambiguous = [], [], []
    for name in names:
        rows = con.execute('SELECT DISTINCT gene_id, kind FROM gene_names WHERE name = ?',
                           (name.upper(),)).fetchall()
        exact = [r['gene_id'] for r in rows if r['kind'] != 'alias']
        ids = exact or sorted({r['gene_id'] for r in rows})
        if not ids:
            prefix, _, mod_id = name.rpartition(':')
            if not prefix or prefix.isalpha():
                ids = [r['id'] for r in con.execute('SELECT id FROM genes WHERE mod_id = ? COLLATE NOCASE',
                                                    (mod_id,))]
        if len(ids) == 1:
            if ids[0] not in resolved:
                resolved.append(ids[0])
        elif ids:
            ambiguous.append({'name': name, 'candidates': list(gene_info(con, ids).values())})
        else:
            unknown.append(name)
    return resolved, unknown, ambiguous


def load_go_index():
    with db() as con:
        if not con.execute("SELECT 1 FROM sqlite_master WHERE name = 'gene_go'").fetchone():
            return None
        info = {r['id']: (r['name'], r['aspect']) for r in con.execute('SELECT * FROM go_terms')}
        gene_terms = {}
        for r in con.execute('SELECT gene_id, term_id FROM gene_go'):
            gene_terms.setdefault(r['gene_id'], set()).add(r['term_id'])
        versions = dict(con.execute('SELECT key, value FROM build_info').fetchall())
        species = {r['id']: r['species'] for r in con.execute('SELECT id, species FROM genes')}
    index = GOIndex(info, gene_terms, species)
    index.versions = {k: v for k, v in versions.items() if k.startswith('go_')}
    return index


def precomputed():
    """Whether the build stores single-gene results (builds made before the setting existed always do)."""
    with db() as con:
        row = con.execute("SELECT value FROM build_info WHERE key = 'precomputed'").fetchone()
    return row is None or row[0] == '1'


def assay_columns():
    """SQL for the assay, value type, display-only and links columns (defaults for builds made before them)."""
    with db() as con:
        cols = {r['name'] for r in con.execute('PRAGMA table_info(datasets)')}
    out = 'assay, value_type' if 'assay' in cols else 'NULL AS assay, NULL AS value_type'
    out += ', display_only, links' if 'display_only' in cols else ', 0 AS display_only, NULL AS links'
    return out


def load_config():
    """The MOD config stored in the build, plus the species present (SGD defaults for older builds)."""
    with db() as con:
        row = con.execute("SELECT value FROM build_info WHERE key = 'config'").fetchone()
        species = [r[0] for r in con.execute(
            'SELECT species FROM datasets WHERE species IS NOT NULL GROUP BY species ORDER BY COUNT(*) DESC')]
    cfg = json.loads(row[0]) if row else modconfig.public(modconfig.DEFAULTS)
    cfg.setdefault('assays', modconfig.DEFAULTS['datasets']['assays'])
    return {**cfg, 'species': species}


@app.get('/api/config')
def config():
    """MOD branding, gene link template, example queries and species, for the front end."""
    return state['config']


@app.get('/api/genes/lookup')
def lookup(q: str = Query(..., description='gene names separated by spaces or commas')):
    names = [n for n in NAME_SEP.split(q) if n][:MAX_QUERY_GENES]
    with db() as con:
        ids, unknown, ambiguous = resolve(con, names)
        info = gene_info(con, ids)
        return {'genes': [info[i] for i in ids], 'unknown': unknown, 'ambiguous': ambiguous}


@app.get('/api/genes/{gene_id}/coexpression')
def coexpression(gene_id: int, limit: int = Query(50, le=100), datasets: int = Query(100, le=100),
                 via: str = Query('', pattern='^(|embed)$')):
    """Single-gene SPELL result: precomputed, or searched live when the build has no precomputed results."""
    if not state['precomputed']:
        return run_live([gene_id], None, limit, datasets, log_kind=via or 'single')
    with db() as con:
        query = gene_info(con, [gene_id]).get(gene_id)
        if query is None:
            raise HTTPException(404, 'unknown gene')
        state['log'].record(via or 'single', [query['systematic']])
        hits = con.execute('SELECT hit_id, score FROM coexpression WHERE gene_id = ? AND rank <= ? '
                           'ORDER BY rank', (gene_id, limit)).fetchall()
        dsets = con.execute('SELECT dataset_id, weight FROM coexpression_datasets WHERE gene_id = ? '
                            'AND rank <= ? ORDER BY rank', (gene_id, datasets)).fetchall()
        ginfo = gene_info(con, [h['hit_id'] for h in hits])
        dinfo = dataset_info(con, [d['dataset_id'] for d in dsets])
    return {
        'query': [query],
        'genes': [{**ginfo[h['hit_id']], 'score': h['score']} for h in hits],
        'datasets': [{**dinfo[d['dataset_id']], 'weight': d['weight']} for d in dsets],
        'source': 'precomputed',
    }


class SearchRequest(BaseModel):
    genes: list[str] = Field(..., min_length=1, max_length=MAX_QUERY_GENES)
    limit: int = Field(50, ge=1, le=100)
    datasets: int = Field(100, ge=1, le=100)
    tags: list[str] = Field(default_factory=list, max_length=50,
                            description='search only datasets carrying any of these tags')
    assay: str | None = Field(None, max_length=50, description='search only datasets of this assay')
    dataset_ids: list[int] | None = Field(None, max_length=10000, description='search only these datasets')


LINK = re.compile(r'\s*(.*?)\s*<(https?://[^>\s]+)>\s*')


def parse_links(value):
    """'GEO GSE1 <https://...> | Atlas <https://...>' -> [{'label', 'url'}] (http(s) links only)."""
    out = []
    for part in (value or '').split(' | '):
        m = LINK.fullmatch(part)
        if m:
            out.append({'label': m.group(1) or m.group(2), 'url': m.group(2)})
    return out


def split_tags(value):
    return [t for t in (value or '').split('|') if t]


def filtered_datasets(con, tags, assay, ids=None):
    """IDs of datasets carrying any of ``tags``, of ``assay`` and among ``ids`` (None: no restriction requested)."""
    if not tags and not assay and ids is None:
        return None
    wanted, allowed = set(tags), None if ids is None else set(ids)
    return [r['id'] for r in con.execute('SELECT id, tags, %s FROM datasets' % state['assay_cols'])
            if (not wanted or wanted & set(split_tags(r['tags']))) and (not assay or r['assay'] == assay)
            and (allowed is None or r['id'] in allowed) and not r['display_only']]


def log_filters(tags, assay, ids=None):
    return list(tags) + (['assay:' + assay] if assay else []) + (['datasets:%d' % len(ids)] if ids is not None else [])


@app.post('/api/search')
def live_search(req: SearchRequest):
    with db() as con:
        ids, unknown, ambiguous = resolve(con, req.genes)
        dataset_ids = filtered_datasets(con, req.tags, req.assay, req.dataset_ids)
    if not ids:
        raise HTTPException(422, {'message': 'no recognised genes', 'unknown': unknown, 'ambiguous': ambiguous})
    if dataset_ids is not None and not dataset_ids:
        raise HTTPException(422, {'message': 'no datasets match the selected tags, assay and datasets'})
    if len(ids) == 1 and dataset_ids is None and state['precomputed']:
        result = coexpression(ids[0], req.limit, req.datasets, '')
        return {**result, 'unknown': unknown, 'ambiguous': ambiguous, 'tags': [], 'assay': None,
                'datasets_searched': None}
    result = run_live(ids, dataset_ids, req.limit, req.datasets, filters=log_filters(req.tags, req.assay,
                                                                                    req.dataset_ids))
    return {**result, 'unknown': unknown, 'ambiguous': ambiguous, 'tags': req.tags, 'assay': req.assay,
            'datasets_searched': len(dataset_ids) if dataset_ids is not None else None}


def run_live(ids, dataset_ids, limit, n_datasets, log_kind=None, filters=()):
    """Search live for SPELL gene IDs, optionally within dataset_ids; returns query, genes and datasets."""
    sem = state['live']
    if not sem.acquire(timeout=30):
        raise HTTPException(503, 'search service busy, try again shortly')
    try:
        result = search(state['comp'], ids, dataset_ids=dataset_ids, max_genes=limit, precision='float64')
    finally:
        sem.release()

    hits = [(g, s) for g, s in zip(result.gene_ids, result.gene_scores) if g not in ids][:limit]
    names = [state['comp'].gene_names[state['comp'].index_of_id[g]] for g in ids]
    state['log'].record(log_kind or ('multi' if len(ids) > 1 else 'single'), names, list(filters))
    weights = result.dataset_weights
    if weights is None:                       # single-gene search: weights only order the datasets
        total = sum(result.raw_weights) or 1.0
        weights = [w / total * 100 for w in result.raw_weights]
    dsets = list(zip(result.dataset_ids, weights))[:n_datasets]
    with db() as con:
        ginfo = gene_info(con, ids + [g for g, _ in hits])
        dinfo = dataset_info(con, [d for d, _ in dsets])
    return {
        'query': [ginfo[i] for i in ids],
        'genes': [{**ginfo[g], 'score': s} for g, s in hits],
        'datasets': [{**dinfo[d], 'weight': w} for d, w in dsets],
        'source': 'live',
    }


@app.get('/api/datasets')
def dataset_list():
    """Every dataset with its publication details, tags and QC flags (for the dataset listing)."""
    with db() as con:
        has_qc = con.execute("SELECT 1 FROM sqlite_master WHERE name = 'dataset_qc'").fetchone()
        qc_cols = ", q.scale, q.warnings" if has_qc else ""
        qc_join = " LEFT JOIN dataset_qc q ON q.dataset_id = d.id" if has_qc else ""
        rows = con.execute('SELECT d.id, d.title, d.first_author, d.authors, d.year, d.journal, d.geo_id, d.pubmed_id, '
                           'd.num_conds, d.num_genes, d.tags, d.has_metadata, d.species, ' +
                           state['assay_cols'] + qc_cols +
                           ' FROM datasets d' + qc_join + ' ORDER BY d.year DESC, d.id DESC').fetchall()
    return [{**dict(r), 'tags': split_tags(r['tags'])} for r in rows]


@app.get('/api/expression-listing')
def expression_listing(genes: str = Query(..., description='gene names'), tags: str = Query(''),
                       assay: str = Query('', description='list only datasets of this assay')):
    """Datasets measuring the query genes, most variable first, without a co-expression search."""
    names = [n for n in NAME_SEP.split(genes) if n][:MAX_QUERY_GENES]
    tag_list = split_tags(tags)
    with db() as con:
        ids, unknown, ambiguous = resolve(con, names)
        if not ids:
            raise HTTPException(422, {'message': 'no recognised genes', 'unknown': unknown, 'ambiguous': ambiguous})
        allowed = filtered_datasets(con, tag_list, assay or None)
        allowed = None if allowed is None else set(allowed)
        query = gene_info(con, ids)
    comp = state['comp']
    rows = []
    for ds in comp.datasets:
        if allowed is not None and ds.id not in allowed:
            continue
        present = [ds.row_of[comp.index_of_id[g]] for g in ids if comp.index_of_id[g] in ds.row_of]
        if not present:
            continue
        with warnings.catch_warnings():         # genes with no values in this dataset give NaN
            warnings.simplefilter('ignore', RuntimeWarning)
            spread = float(np.nanmean(np.nanvar(ds.values[present], axis=1)))
        rows.append((ds.id, len(present), 0.0 if np.isnan(spread) else spread))
    rows.sort(key=lambda r: (-r[1], -r[2]))
    with db() as con:
        dinfo = dataset_info(con, [r[0] for r in rows])
    state['log'].record('listing', [query[i]['systematic'] for i in ids], log_filters(tag_list, assay or None))
    return {
        'query': [query[i] for i in ids], 'unknown': unknown, 'ambiguous': ambiguous, 'tags': tag_list,
        'assay': assay or None,
        'datasets': [{**dinfo[d], 'genes_measured': n, 'variance': v} for d, n, v in rows],
    }


@app.get('/api/datasets/{dataset_id}/download')
def dataset_download(dataset_id: int):
    """The dataset's expression matrix as searched (after QC), as tab-separated text."""
    ds = next((d for d in state['comp'].datasets if d.id == dataset_id), None)
    if ds is None:
        raise HTTPException(404, 'unknown dataset')
    comp = state['comp']
    with db() as con:
        std = {r['id']: r['standard'] for r in con.execute('SELECT id, standard FROM genes')}
        geo = con.execute('SELECT geo_id FROM datasets WHERE id = ?', (dataset_id,)).fetchone()['geo_id']
    def clean(text):
        return str(text).replace('\t', ' ').replace('\n', ' ')

    lines = ['\t'.join(['systematic_name', 'standard_name'] + [clean(c) for c in ds.conditions])]
    for r, gi in enumerate(ds.gene_index):
        gid = comp.gene_ids[gi]
        vals = ['' if np.isnan(v) else '%.4g' % v for v in ds.values[r]]
        lines.append('\t'.join([comp.gene_names[gi], std.get(gid) or ''] + vals))
    name = '%s_spell-next_%d.tsv' % (geo or 'dataset', dataset_id)
    return PlainTextResponse('\n'.join(lines) + '\n', media_type='text/tab-separated-values',
                             headers={'Content-Disposition': 'attachment; filename="%s"' % name})


@app.get('/api/info')
def info():
    """Counts and data versions for the help and about pages."""
    with db() as con:
        has = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}
        versions = {k: v for k, v in con.execute('SELECT key, value FROM build_info') if k.startswith('go_')} \
            if 'build_info' in has else {}
        qc = dict(con.execute("SELECT scale, COUNT(*) FROM dataset_qc GROUP BY scale").fetchall()) \
            if 'dataset_qc' in has else {}
        flagged = con.execute("SELECT COUNT(*) FROM dataset_qc WHERE warnings != ''").fetchone()[0] \
            if 'dataset_qc' in has else 0
        searched = ('(SELECT id, has_metadata, num_conds, %s FROM datasets) WHERE display_only = 0'
                    % state['assay_cols'])
        return {
            'datasets': con.execute('SELECT COUNT(*) FROM ' + searched).fetchone()[0],
            'datasets_with_metadata': con.execute('SELECT SUM(has_metadata) FROM ' + searched).fetchone()[0],
            'datasets_display_only': con.execute('SELECT COUNT(*) FROM (SELECT %s FROM datasets) WHERE display_only = 1'
                                                 % state['assay_cols']).fetchone()[0],
            'conditions': con.execute('SELECT SUM(num_conds) FROM ' + searched).fetchone()[0],
            'geo_series': con.execute("SELECT COUNT(DISTINCT geo_id) FROM datasets WHERE geo_id != '' AND id IN "
                                      '(SELECT id FROM ' + searched + ')').fetchone()[0],
            'publications': con.execute('SELECT COUNT(DISTINCT pubmed_id) FROM datasets WHERE pubmed_id AND id IN '
                                        '(SELECT id FROM ' + searched + ')').fetchone()[0],
            'genes': con.execute('SELECT COUNT(*) FROM genes').fetchone()[0],
            'qc_scales': qc,
            'qc_flagged': flagged,
            'versions': versions,
            'max_query_genes': MAX_QUERY_GENES,
            'precomputed': state['precomputed'],
        }


DOWNLOADS = {
    'coexpression_genes.tsv.gz': (
        ['query_systematic', 'query_name', 'rank', 'gene_systematic', 'gene_name', 'score'],
        'SELECT q.systematic, q.standard, c.rank, h.systematic, h.standard, c.score FROM coexpression c '
        'JOIN genes q ON q.id = c.gene_id JOIN genes h ON h.id = c.hit_id ORDER BY q.systematic, c.rank'),
    'coexpression_datasets.tsv.gz': (
        ['query_systematic', 'query_name', 'rank', 'dataset_id', 'geo_id', 'weight_percent', 'dataset_title'],
        'SELECT q.systematic, q.standard, c.rank, d.id, d.geo_id, c.weight, d.title FROM coexpression_datasets c '
        'JOIN genes q ON q.id = c.gene_id JOIN datasets d ON d.id = c.dataset_id ORDER BY q.systematic, c.rank'),
}


@app.get('/api/downloads/{name}')
def download(name: str):
    """Every gene's precomputed top co-expressed genes or datasets, as gzipped tab-separated text."""
    if name not in DOWNLOADS:
        raise HTTPException(404, 'unknown download')
    if not state['precomputed']:
        raise HTTPException(404, 'this site has no precomputed results (searches run live)')
    cache = state.setdefault('downloads', {})
    if name not in cache:
        cols, sql = DOWNLOADS[name]
        buf = io.StringIO()
        buf.write('\t'.join(cols) + '\n')
        with db() as con:
            for row in con.execute(sql):
                buf.write('\t'.join('%.4g' % v if isinstance(v, float) else
                                    str(v if v is not None else '').replace('\t', ' ') for v in row) + '\n')
        cache[name] = gzip.compress(buf.getvalue().encode(), 6)
    return Response(cache[name], media_type='application/gzip',
                    headers={'Content-Disposition': 'attachment; filename="%s_%s"' % (
                        state['config']['mod']['id'].lower(), name)})


@app.get('/api/assays')
def assay_list():
    """Assays with the number of datasets of each (datasets with no assay recorded are not counted)."""
    with db() as con:
        rows = con.execute('SELECT assay, COUNT(*) AS n FROM (SELECT %s FROM datasets) WHERE assay IS NOT NULL '
                           'GROUP BY assay ORDER BY n DESC, assay' % state['assay_cols']).fetchall()
    return [{'assay': r['assay'], 'count': r['n']} for r in rows]


@app.get('/api/tags')
def tag_list():
    """Dataset tags with the number of datasets carrying each, and their keyword category, definition and
    link (null for tags not in the keyword list, and for builds without one)."""
    counts = {}
    with db() as con:
        for (tags,) in con.execute('SELECT tags FROM datasets'):
            for t in split_tags(tags):
                counts[t] = counts.get(t, 0) + 1
        has_kw = con.execute("SELECT 1 FROM sqlite_master WHERE name = 'keywords'").fetchone()
        kw = {r['name']: dict(r) for r in con.execute('SELECT * FROM keywords')} if has_kw else {}
    return [{'tag': t, 'count': n, 'category': kw.get(t, {}).get('category'),
             'definition': kw.get(t, {}).get('definition'), 'link': kw.get(t, {}).get('link'),
             'position': kw.get(t, {}).get('position')}
            for t, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))]


class EnrichmentRequest(BaseModel):
    gene_ids: list[int] = Field(..., min_length=1, max_length=500)
    aspect: str = Field('P', pattern='^[PFC]$')
    max_results: int = Field(25, ge=1, le=200)


@app.post('/api/go-enrichment')
def go_enrichment(req: EnrichmentRequest):
    """GO terms over-represented among the given genes, against all annotated genes."""
    index = state.get('go')
    if index is None:
        raise HTTPException(503, 'GO annotations not loaded')
    results = index.enrich(req.gene_ids, req.aspect, max_results=req.max_results)
    species = index.species_of(req.gene_ids)
    with db() as con:
        names = {g['id']: g['standard'] or g['systematic']
                 for g in gene_info(con, sorted({g for r in results for g in r.genes})).values()}
    return {
        'aspect': ASPECTS[req.aspect],
        'species': species,
        'list_size': results[0].list_size if results else 0,
        'background_size': len(index.backgrounds[species][req.aspect]) if species in index.backgrounds else 0,
        'versions': index.versions,
        'terms': [{**r.__dict__, 'genes': [names.get(g, str(g)) for g in r.genes],
                   'also': [{'term_id': t, 'name': n} for t, n in r.also]} for r in results],
    }


@app.get('/api/datasets/{dataset_id}')
def dataset_detail(dataset_id: int):
    with db() as con:
        row = con.execute('SELECT * FROM datasets WHERE id = ?', (dataset_id,)).fetchone()
        if row is None:
            raise HTTPException(404, 'unknown dataset')
        has_qc = con.execute("SELECT 1 FROM sqlite_master WHERE name = 'dataset_qc'").fetchone()
        qc = has_qc and con.execute('SELECT scale, actions, warnings, cond_median_spread, z_mean_before, '
                                    'z_mean_after FROM dataset_qc WHERE dataset_id = ?', (dataset_id,)).fetchone()
    ds = next((d for d in state['comp'].datasets if d.id == dataset_id), None)
    return {**dict(row), 'display_only': int(ds is None), 'links': parse_links(row['links'] if 'links' in row.keys()
                                                                               else None),
            'conditions': ds.conditions if ds else [], 'qc': dict(qc) if qc and ds else None}


@app.get('/api/datasets/{dataset_id}/expression')
def dataset_expression(dataset_id: int, genes: str = Query(..., description='comma-separated SPELL gene IDs'),
                       raw: bool = Query(False, description='also return values as stored (after QC)')):
    """Expression values for the requested genes: centred on each gene's mean (for heatmaps)
    and, with ``raw=true``, the stored values too."""
    ds = next((d for d in state['comp'].datasets if d.id == dataset_id), None)
    if ds is None:
        raise HTTPException(404, 'unknown dataset')
    try:
        wanted = [int(g) for g in genes.split(',') if g][:200]
    except ValueError:
        raise HTTPException(422, 'genes must be numeric IDs')
    comp = state['comp']

    def as_list(v):
        return [None if np.isnan(x) else round(float(x), 3) for x in v]

    rows = []
    for gid in wanted:
        row = ds.row_of.get(comp.index_of_id.get(gid, -1))
        if row is None:
            rows.append({'gene_id': gid, 'values': None, **({'raw': None} if raw else {})})
            continue
        v = ds.values[row]
        centred = v - np.nanmean(v) if np.any(~np.isnan(v)) else v
        rows.append({'gene_id': gid, 'values': as_list(centred), **({'raw': as_list(v)} if raw else {})})
    return {'dataset_id': dataset_id, 'conditions': ds.conditions, 'rows': rows}


app.mount('/static', StaticFiles(directory=STATIC_DIR), name='static')


@app.get('/embed')
def embed():
    """Compact single-gene result for other sites to show in an iframe: /embed?gene=ACT1&genes=10&sets=5"""
    return FileResponse(STATIC_DIR / 'embed.html')


@app.get('/')
def index():
    return FileResponse(STATIC_DIR / 'index.html')
