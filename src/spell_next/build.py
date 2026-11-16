"""Build the spell-next data store for one MOD from a copy of its SPELL files.

Usage: python -m spell_next.build RAW_DIR OUT_DIR [--workers N] [--limit N]

RAW_DIR holds gene_list.txt, dataset_list.txt, expression/, an optional mod.toml
(see spell_next.modconfig; SGD defaults otherwise) and the gene, dataset metadata
and GO files it names. OUT_DIR receives:

* compendium.npz - expression matrices, parsed correctly (no tab collapsing,
  duplicate dataset entries dropped) and normalized by spell_next.qc;
* spell.sqlite  - gene and dataset metadata, QC records and precomputed single-gene results;
* qc_report.tsv - per-dataset QC summary for curators.
"""
import argparse
import csv
import gzip
import io
import json
import os
import re
import sqlite3
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from spell_next import go, modconfig, qc
from spell_next.compendium import Compendium
from spell_next.search import search

TOP_GENES = 100
TOP_DATASETS = 100

SCHEMA = """
CREATE TABLE genes (
    id INTEGER PRIMARY KEY,         -- SPELL gene ID
    systematic TEXT NOT NULL,
    standard TEXT,
    mod_id TEXT,                    -- the MOD's gene ID (SGD S000001855, WormBase WBGene00006789)
    qualifier TEXT,
    description TEXT,
    species TEXT
);
CREATE TABLE gene_names (
    name TEXT NOT NULL,             -- upper case systematic, standard or alias name
    gene_id INTEGER NOT NULL REFERENCES genes(id),
    kind TEXT NOT NULL              -- systematic | standard | alias
);
CREATE INDEX gene_names_name ON gene_names(name);
CREATE TABLE datasets (
    id INTEGER PRIMARY KEY,         -- SPELL dataset ID
    filename TEXT NOT NULL,
    geo_id TEXT,
    pubmed_id INTEGER,
    title TEXT,
    description TEXT,
    first_author TEXT,
    authors TEXT,
    journal TEXT,
    year INTEGER,
    tags TEXT,
    num_conds INTEGER,
    num_genes INTEGER,
    has_metadata INTEGER NOT NULL,
    species TEXT,                   -- species of most of the dataset's genes
    assay TEXT,                     -- e.g. microarray | RNA-seq (datasets.assays in mod.toml)
    value_type TEXT                 -- what the values are, e.g. log2 TMM CPM, log2 TPM, log2 ratio
);
CREATE TABLE dataset_qc (
    dataset_id INTEGER PRIMARY KEY REFERENCES datasets(id),
    scale TEXT,                     -- log | ratio | intensity | unclear | empty
    actions TEXT,                   -- what normalization did, '; '-separated
    warnings TEXT,                  -- issues for a curator to review, '; '-separated
    cond_median_spread REAL,        -- before normalization (log units)
    z_mean_before REAL,             -- mean Fisher z over all gene pairs, before ...
    z_mean_after REAL,              -- ... and after normalization
    z_std_after REAL
);
CREATE TABLE go_terms (
    id TEXT PRIMARY KEY,            -- GO:0006096
    name TEXT NOT NULL,
    aspect TEXT NOT NULL            -- P | F | C
);
CREATE TABLE gene_go (
    gene_id INTEGER NOT NULL,
    term_id TEXT NOT NULL,          -- annotation propagated to all ancestors (is_a, part_of)
    PRIMARY KEY (gene_id, term_id)
) WITHOUT ROWID;
CREATE TABLE keywords (            -- the MOD's tag vocabulary, in the order of its keyword list
    name TEXT PRIMARY KEY,
    category TEXT,                  -- groups tags in the tag picker
    definition TEXT,
    link TEXT,                      -- the keyword's page on the MOD's site
    position INTEGER NOT NULL
);
CREATE TABLE build_info (
    key TEXT PRIMARY KEY,
    value TEXT
);
CREATE TABLE coexpression (
    gene_id INTEGER NOT NULL,
    rank INTEGER NOT NULL,          -- 1 = best hit (the query gene itself is excluded)
    hit_id INTEGER NOT NULL,
    score REAL NOT NULL,            -- sqrt of the SPELL score, as SPELL displays it
    PRIMARY KEY (gene_id, rank)
) WITHOUT ROWID;
CREATE TABLE coexpression_datasets (
    gene_id INTEGER NOT NULL,
    rank INTEGER NOT NULL,          -- 1 = most informative dataset for this gene
    dataset_id INTEGER NOT NULL,
    weight REAL NOT NULL,           -- percent of total weight
    PRIMARY KEY (gene_id, rank)
) WITHOUT ROWID;
"""


def read_sgd_features(path, default_species=None):
    """Gene info from SGD's SGD_features.tab: {systematic name: gene dict}."""
    out = {}
    with open(path, errors='replace') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 16 or not f[3] or f[1] in ('CDS', 'intron', 'five_prime_UTR_intron'):
                continue
            out.setdefault(f[3], {'mod_id': f[0], 'standard': f[4] or None,
                                  'aliases': [a for a in f[5].split('|') if a],
                                  'qualifier': f[2] or None, 'description': f[15] or None,
                                  'species': default_species})
    return out


def read_gene_table(path, default_species=None):
    """Gene info from a generic tab-separated table with a header row.

    Columns (only ``systematic`` is required): systematic, mod_id, standard,
    aliases ('|'-separated), species, description, qualifier.
    """
    out = {}
    with open(path, newline='', errors='replace') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            name = (row.get('systematic') or '').strip()
            if not name:
                continue
            out.setdefault(name, {
                'mod_id': (row.get('mod_id') or '').strip() or None,
                'standard': (row.get('standard') or '').strip() or None,
                'aliases': [a.strip() for a in (row.get('aliases') or '').split('|') if a.strip()],
                'qualifier': (row.get('qualifier') or '').strip() or None,
                'description': (row.get('description') or '').strip() or None,
                'species': (row.get('species') or '').strip() or default_species,
            })
    return out


def read_dataset_dump(path):
    """Dataset metadata from SPELL's SQLdump*_fordownloadfiles.txt export: {filename: dict}."""
    out = {}
    with open(path, errors='replace') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 11:
                continue
            out[f[1]] = {
                'pubmed_id': int(f[0]) if f[0].isdigit() else None,
                'title': f[2] or f[7], 'description': f[3],
                'first_author': f[5], 'authors': f[6], 'journal': f[8],
                'year': int(f[9]) if f[9].isdigit() else None,
                'tags': f[10], 'geo_id': None,
            }
    return out


def read_dataset_table(path):
    """Dataset metadata from a tab-separated table with a header row: {filename: dict}.

    Columns (only ``filename`` is required): filename, geo_id, pubmed_id, title,
    description, first_author, authors, journal, year, tags ('|'-separated), assay,
    value_type.
    """
    out = {}
    with open(path, newline='', errors='replace') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            name = (row.get('filename') or '').strip()
            if not name:
                continue
            def get(key, row=row):
                return (row.get(key) or '').strip() or None
            out[name] = {
                'pubmed_id': int(get('pubmed_id')) if (get('pubmed_id') or '').isdigit() else None,
                'year': int(get('year')) if (get('year') or '').isdigit() else None,
                'geo_id': get('geo_id'), 'title': get('title'), 'description': get('description'),
                'first_author': get('first_author'), 'authors': get('authors'),
                'journal': get('journal'), 'tags': get('tags'),
                'assay': get('assay'), 'value_type': get('value_type'),
            }
    return out


def read_keywords(path):
    """The keyword list (comma- or tab-separated, header row) as a list of dicts, in file order.

    Columns: display_name (or name/tag/keyword; required), category, description (or
    definition), format_name (used in keyword links).
    """
    with open(path, newline='', encoding='utf-8-sig', errors='replace') as fh:
        text = fh.read()
    delimiter = '\t' if '\t' in text.split('\n', 1)[0] else ','
    out, seen = [], set()
    for row in csv.DictReader(io.StringIO(text), delimiter=delimiter):
        row = {(k or '').strip().lower(): (v or '').strip() for k, v in row.items()}
        name = row.get('display_name') or row.get('name') or row.get('tag') or row.get('keyword')
        if not name or name in seen:
            continue
        seen.add(name)
        out.append({'name': name, 'category': row.get('category') or None,
                    'definition': row.get('description') or row.get('definition') or None,
                    'format_name': row.get('format_name') or name.replace(' ', '_')})
    return out


def write_keywords(db, raw_dir, cfg):
    rel = cfg['datasets'].get('keywords')
    path = modconfig.path(raw_dir, rel)
    if not path:
        return
    if not os.path.exists(path):
        print('keywords: %s not found; tags are shown without categories' % rel)
        return
    keywords = read_keywords(path)
    template = cfg['links'].get('keyword')
    db.executemany('INSERT INTO keywords VALUES (?,?,?,?,?)', [
        (k['name'], k['category'], k['definition'],
         template.format(format_name=k['format_name'], name=k['name']) if template else None, i)
        for i, k in enumerate(keywords)])
    print('keywords: %d in %d categories' % (len(keywords), len({k['category'] for k in keywords if k['category']})))


def read_gene_info(raw_dir, cfg):
    g = cfg['genes']
    path = modconfig.path(raw_dir, g['file'])
    if not path or not os.path.exists(path):
        print('genes: %s not found; using systematic names only' % g['file'])
        return {}
    reader = read_sgd_features if g['format'] == 'sgd_features' else read_gene_table
    return reader(path, g.get('default_species'))


def read_dataset_info(raw_dir, cfg):
    d = cfg['datasets']
    path = modconfig.path(raw_dir, d['metadata'])
    if not path or not os.path.exists(path):
        print('datasets: no metadata file; datasets get titles from file names')
        return {}
    return read_dataset_dump(path) if d['format'] == 'spell_dump' else read_dataset_table(path)


def write_metadata(db, comp, genes, datasets, default_species=None):
    db.executescript(SCHEMA)
    gene_rows, name_rows, species_of = [], [], {}
    for gid, sys_name in zip(comp.gene_ids, comp.gene_names):
        info = genes.get(sys_name, {})
        species = info.get('species') or default_species
        species_of[gid] = species
        gene_rows.append((gid, sys_name, info.get('standard'), info.get('mod_id'), info.get('qualifier'),
                          info.get('description'), species))
        name_rows.append((sys_name.upper(), gid, 'systematic'))
        if info.get('standard'):
            name_rows.append((info['standard'].upper(), gid, 'standard'))
        name_rows.extend((a.upper(), gid, 'alias') for a in info.get('aliases', []))
    db.executemany('INSERT INTO genes VALUES (?,?,?,?,?,?,?)', gene_rows)
    db.executemany('INSERT INTO gene_names VALUES (?,?,?)', name_rows)

    for ds in comp.datasets:
        meta = datasets.get(ds.filename)
        geo = (meta and meta.get('geo_id')) or (re.match(r'(GSE\d+)', ds.filename) or [None, None])[1]
        counts = Counter(species_of[comp.gene_ids[i]] for i in ds.gene_index)
        species = counts.most_common(1)[0][0] if counts else default_species
        db.execute('INSERT INTO datasets VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', (
            ds.id, ds.filename, geo, meta and meta['pubmed_id'],
            (meta and meta['title']) or geo or ds.filename,
            meta and meta['description'], meta and meta['first_author'], meta and meta['authors'],
            meta and meta['journal'], meta and meta['year'], meta and meta['tags'],
            ds.num_conds, len(ds.gene_index), int(bool(meta and meta.get('title'))), species,
            meta and meta.get('assay'), meta and meta.get('value_type')))


def _qc_one(values, centre=False):
    diag = qc.diagnose(values)
    scale = 'empty' if diag['empty'] else qc.classify_scale(diag)
    new_values, actions, warnings = qc.normalize(values, centre=centre)
    # Report condition shifts on the corrected scale (log units), before any centring.
    rescaled = new_values if not centre else qc.normalize(values)[0]
    spread = qc.diagnose(rescaled).get('cond_median_spread')
    z_after, z_std = qc.pair_stats(new_values)
    # Most datasets are left unchanged; then the "before" statistics are the same (and costly).
    unchanged = np.array_equal(new_values, values, equal_nan=True)
    z_before = z_after if unchanged else qc.pair_stats(values)[0]
    return new_values, scale, actions, warnings, spread, z_before, z_after, z_std


def apply_qc(db, comp, workers, centre=False):
    """Normalize every dataset in place, recompute its z statistics and record what was done."""
    t0 = time.time()
    with ProcessPoolExecutor(workers) as pool:
        results = pool.map(_qc_one, [d.values for d in comp.datasets], [centre] * len(comp.datasets),
                           chunksize=4)
        for ds, (values, scale, actions, warnings, spread, z_before, z_after, z_std) in zip(comp.datasets, results):
            ds.values = values
            ds.z_mean, ds.z_std = z_after, z_std
            db.execute('INSERT INTO dataset_qc VALUES (?,?,?,?,?,?,?,?)',
                       (ds.id, scale, '; '.join(actions), '; '.join(warnings), spread, z_before, z_after, z_std))
    db.commit()
    print('qc: normalized %d datasets (%.0fs)' % (len(comp.datasets), time.time() - t0))


def write_qc_report(db, path):
    rows = db.execute('SELECT d.id, d.filename, d.title, q.scale, q.cond_median_spread, q.z_mean_before, '
                      'q.z_mean_after, q.actions, q.warnings FROM dataset_qc q JOIN datasets d ON d.id = q.dataset_id '
                      'ORDER BY (q.warnings != \'\') DESC, q.z_mean_before - q.z_mean_after DESC').fetchall()
    with open(path, 'w') as out:
        out.write('dataset_id\tfilename\ttitle\tscale\tcond_median_spread\tz_mean_before\tz_mean_after\t'
                  'actions\twarnings\n')
        for r in rows:
            out.write('\t'.join('' if v is None else ('%.3f' % v if isinstance(v, float) else str(v))
                                for v in r) + '\n')


def write_go_tables(db, raw_dir, cfg):
    """Store GO terms and propagated annotations (matched on the MOD gene ID) for the compendium's genes."""
    obo = modconfig.path(raw_dir, cfg['go']['obo'])
    gafs = [modconfig.path(raw_dir, g) for g in cfg['go']['gaf']]
    missing = [p for p in [obo] + gafs if not p or not os.path.exists(p)]
    if missing:
        print('go: %s missing, GO enrichment disabled' % ', '.join(missing))
        return
    terms = go.read_obo(obo)
    mod_to_gene = {m: g for g, m in db.execute('SELECT id, mod_id FROM genes WHERE mod_id IS NOT NULL')}
    direct = {}
    for gaf in gafs:
        for gene, ts in go.read_gaf(gaf, mod_to_gene, terms).items():
            direct.setdefault(gene, set()).update(ts)
    gene_terms = go.propagate(direct, go.ancestors(terms))
    used = set().union(*gene_terms.values()) if gene_terms else set()
    db.executemany('INSERT INTO go_terms VALUES (?,?,?)', [(t, terms[t][0], terms[t][1]) for t in sorted(used)])
    db.executemany('INSERT INTO gene_go VALUES (?,?)', [(g, t) for g, ts in gene_terms.items() for t in ts])
    with open(obo) as fh:
        version = next((ln.split(': ', 1)[1].strip() for ln in fh if ln.startswith('data-version')), 'unknown')
    dates = []
    for gaf in gafs:
        opener = gzip.open if gaf.endswith('.gz') else open
        with opener(gaf, 'rt') as fh:
            dates.append(next((ln.split(':', 1)[1].strip() for ln in fh if ln.startswith('!date-generated')
                               or ln.startswith('!Generated')), 'unknown'))
    db.executemany('INSERT OR REPLACE INTO build_info VALUES (?,?)',
                   [('go_ontology', version), ('go_annotations', ', '.join(dates))])
    db.commit()
    print('go: %d genes annotated, %d terms (ontology %s)' % (len(gene_terms), len(used), version))


_worker_comp = None


def _init_worker(npz_path):
    global _worker_comp
    _worker_comp = Compendium.load(npz_path)


def _search_one(gene_id):
    result = search(_worker_comp, [gene_id], max_genes=TOP_GENES, precision='float64')
    hits = [(g, s) for g, s in zip(result.gene_ids, result.gene_scores) if g != gene_id][:TOP_GENES]
    total = sum(result.raw_weights) or 1.0
    dsets = [(d, w / total * 100) for d, w in zip(result.dataset_ids, result.raw_weights)][:TOP_DATASETS]
    return gene_id, hits, dsets


def precompute(db, comp, npz_path, gene_ids, workers):
    """Store the top genes and datasets for every gene, computed dataset by dataset."""
    from spell_next.precompute import precompute_all

    t0 = time.time()
    results = precompute_all(comp, npz_path, TOP_GENES, TOP_DATASETS, workers)
    for gene_id in gene_ids:
        hits, dsets = results[gene_id]
        db.executemany('INSERT INTO coexpression VALUES (?,?,?,?)',
                       [(gene_id, r, g, s) for r, (g, s) in enumerate(hits, 1)])
        db.executemany('INSERT INTO coexpression_datasets VALUES (?,?,?,?)',
                       [(gene_id, r, d, w) for r, (d, w) in enumerate(dsets, 1)])
    db.commit()
    print('  stored %d genes (%.0fs)' % (len(gene_ids), time.time() - t0))


def precompute_per_gene(db, npz_path, gene_ids, workers):
    """Reference implementation: one search per gene (slow; used to check precompute)."""
    t0, done = time.time(), 0
    with ProcessPoolExecutor(workers, initializer=_init_worker, initargs=(npz_path,)) as pool:
        for gene_id, hits, dsets in pool.map(_search_one, gene_ids, chunksize=8):
            db.executemany('INSERT INTO coexpression VALUES (?,?,?,?)',
                           [(gene_id, r, g, s) for r, (g, s) in enumerate(hits, 1)])
            db.executemany('INSERT INTO coexpression_datasets VALUES (?,?,?,?)',
                           [(gene_id, r, d, w) for r, (d, w) in enumerate(dsets, 1)])
            done += 1
            if done % 500 == 0:
                db.commit()
                print('  %d/%d genes, %.0fs' % (done, len(gene_ids), time.time() - t0), flush=True)
    db.commit()


def run_build(raw_dir, out_dir, workers=None, limit=None, no_qc=False, centre=False, per_gene=False):
    """Build OUT_DIR from RAW_DIR; returns the path of the SQLite file."""
    workers = workers or os.cpu_count()
    os.makedirs(out_dir, exist_ok=True)
    npz_path = os.path.join(out_dir, 'compendium.npz')
    db_path = os.path.join(out_dir, 'spell.sqlite')

    t0 = time.time()
    cfg = modconfig.load(raw_dir)
    print('mod: %s (%s)' % (cfg['mod']['id'], cfg['mod']['name']))
    comp = Compendium.from_spell_files(raw_dir, faithful=False, dedupe=True)
    print('compendium: %d datasets, %d genes (%.0fs)' % (len(comp.datasets), len(comp.gene_ids), time.time() - t0))

    if os.path.exists(db_path):
        os.remove(db_path)
    db = sqlite3.connect(db_path)
    write_metadata(db, comp, read_gene_info(raw_dir, cfg), read_dataset_info(raw_dir, cfg),
                   cfg['genes'].get('default_species'))
    write_keywords(db, raw_dir, cfg)
    db.executemany('INSERT OR REPLACE INTO build_info VALUES (?,?)',
                   [('config', json.dumps(modconfig.public(cfg))), ('mod', cfg['mod']['id']),
                    ('built', time.strftime('%Y-%m-%d %H:%M:%S'))])
    db.commit()
    n_meta = db.execute('SELECT SUM(has_metadata), COUNT(*) FROM datasets').fetchone()
    print('metadata: %d/%d datasets have publication metadata' % n_meta)
    species = db.execute('SELECT species, COUNT(*) FROM datasets GROUP BY species ORDER BY 2 DESC').fetchall()
    print('species: ' + ', '.join('%s (%d datasets)' % s for s in species))

    write_go_tables(db, raw_dir, cfg)

    if not no_qc:
        apply_qc(db, comp, workers, centre=centre)
        write_qc_report(db, os.path.join(out_dir, 'qc_report.tsv'))
    else:
        # Datasets without a .sinfo file still need their z statistics.
        for ds in comp.datasets:
            if ds.z_mean != ds.z_mean or ds.z_std != ds.z_std:
                ds.z_mean, ds.z_std = qc.pair_stats(ds.values)
    comp.save(npz_path)

    gene_ids = comp.gene_ids[:limit] if limit else comp.gene_ids
    print('precomputing %d genes with %d workers...' % (len(gene_ids), workers))
    if per_gene:
        precompute_per_gene(db, npz_path, gene_ids, workers)
    else:
        precompute(db, comp, npz_path, gene_ids, workers)
    db.execute('VACUUM')
    db.close()
    print('done in %.0fs: %s' % (time.time() - t0, db_path))
    return db_path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('raw_dir')
    ap.add_argument('out_dir')
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--limit', type=int, help='precompute only the first N genes (for testing)')
    ap.add_argument('--no-qc', action='store_true', help='skip normalization (use values as published)')
    ap.add_argument('--centre', action='store_true',
                    help='also centre each condition on its median (off by default: lowers the GO benchmark)')
    ap.add_argument('--per-gene', action='store_true',
                    help='precompute with one search per gene (slow reference implementation)')
    args = ap.parse_args()
    run_build(args.raw_dir, args.out_dir, args.workers, args.limit, args.no_qc, args.centre, args.per_gene)


if __name__ == '__main__':
    main()
