"""spell-next web API and front end.

Run: SPELL_NEXT_DATA=data/build uvicorn spell_next.api:app

Single-gene results are served from the precomputed table; multi-gene searches run
live against the in-memory compendium, a few at a time.
"""
import os
import re
import sqlite3
import threading
from contextlib import asynccontextmanager
from pathlib import Path

import numpy as np
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from spell_next.compendium import Compendium
from spell_next.go import ASPECTS, GOIndex
from spell_next.search import search

DATA_DIR = Path(os.environ.get('SPELL_NEXT_DATA', 'data/build'))
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
    rows = con.execute('SELECT id, geo_id, pubmed_id, title, first_author, year, tags, num_conds, num_genes '
                       'FROM datasets WHERE id IN (%s)' % ','.join('?' * len(ids)), list(ids))
    return {r['id']: dict(r) for r in rows}


def resolve(con, names):
    """Resolve gene names to SPELL gene IDs: systematic and standard names first, then aliases."""
    resolved, unknown, ambiguous = [], [], []
    for name in names:
        rows = con.execute('SELECT DISTINCT gene_id, kind FROM gene_names WHERE name = ?',
                           (name.upper(),)).fetchall()
        exact = [r['gene_id'] for r in rows if r['kind'] != 'alias']
        ids = exact or sorted({r['gene_id'] for r in rows})
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
    index = GOIndex(info, gene_terms)
    index.versions = versions
    return index


@app.get('/api/genes/lookup')
def lookup(q: str = Query(..., description='gene names separated by spaces or commas')):
    names = [n for n in NAME_SEP.split(q) if n][:MAX_QUERY_GENES]
    with db() as con:
        ids, unknown, ambiguous = resolve(con, names)
        info = gene_info(con, ids)
        return {'genes': [info[i] for i in ids], 'unknown': unknown, 'ambiguous': ambiguous}


@app.get('/api/genes/{gene_id}/coexpression')
def coexpression(gene_id: int, limit: int = Query(50, le=100), datasets: int = Query(20, le=100)):
    """Precomputed single-gene SPELL result."""
    with db() as con:
        query = gene_info(con, [gene_id]).get(gene_id)
        if query is None:
            raise HTTPException(404, 'unknown gene')
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
    datasets: int = Field(20, ge=1, le=100)


@app.post('/api/search')
def live_search(req: SearchRequest):
    with db() as con:
        ids, unknown, ambiguous = resolve(con, req.genes)
    if not ids:
        raise HTTPException(422, {'message': 'no recognised genes', 'unknown': unknown, 'ambiguous': ambiguous})
    if len(ids) == 1:
        result = coexpression(ids[0], req.limit, req.datasets)
        return {**result, 'unknown': unknown, 'ambiguous': ambiguous}

    sem = state['live']
    if not sem.acquire(timeout=30):
        raise HTTPException(503, 'search service busy, try again shortly')
    try:
        result = search(state['comp'], ids, max_genes=req.limit, precision='float64')
    finally:
        sem.release()

    hits = [(g, s) for g, s in zip(result.gene_ids, result.gene_scores) if g not in ids][:req.limit]
    dsets = list(zip(result.dataset_ids, result.dataset_weights))[:req.datasets]
    with db() as con:
        ginfo = gene_info(con, ids + [g for g, _ in hits])
        dinfo = dataset_info(con, [d for d, _ in dsets])
    return {
        'query': [ginfo[i] for i in ids],
        'genes': [{**ginfo[g], 'score': s} for g, s in hits],
        'datasets': [{**dinfo[d], 'weight': w} for d, w in dsets],
        'unknown': unknown, 'ambiguous': ambiguous, 'source': 'live',
    }


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
    with db() as con:
        names = {g['id']: g['standard'] or g['systematic']
                 for g in gene_info(con, sorted({g for r in results for g in r.genes})).values()}
    return {
        'aspect': ASPECTS[req.aspect],
        'list_size': results[0].list_size if results else 0,
        'background_size': len(index.background[req.aspect]),
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
    ds = next(d for d in state['comp'].datasets if d.id == dataset_id)
    return {**dict(row), 'conditions': ds.conditions, 'qc': dict(qc) if qc else None}


@app.get('/api/datasets/{dataset_id}/expression')
def dataset_expression(dataset_id: int, genes: str = Query(..., description='comma-separated SPELL gene IDs')):
    """Expression values (row-centred) for the requested genes, for heatmaps."""
    ds = next((d for d in state['comp'].datasets if d.id == dataset_id), None)
    if ds is None:
        raise HTTPException(404, 'unknown dataset')
    try:
        wanted = [int(g) for g in genes.split(',') if g][:200]
    except ValueError:
        raise HTTPException(422, 'genes must be numeric IDs')
    comp = state['comp']
    rows = []
    for gid in wanted:
        row = ds.row_of.get(comp.index_of_id.get(gid, -1))
        if row is None:
            rows.append({'gene_id': gid, 'values': None})
            continue
        v = ds.values[row]
        centred = v - np.nanmean(v) if np.any(~np.isnan(v)) else v
        rows.append({'gene_id': gid, 'values': [None if np.isnan(x) else round(float(x), 3) for x in centred]})
    return {'dataset_id': dataset_id, 'conditions': ds.conditions, 'rows': rows}


app.mount('/static', StaticFiles(directory=STATIC_DIR), name='static')


@app.get('/')
def index():
    return FileResponse(STATIC_DIR / 'index.html')
