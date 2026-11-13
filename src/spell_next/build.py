"""Build the spell-next data store from a copy of SPELL's files.

Usage: python -m spell_next.build RAW_DIR OUT_DIR [--workers N] [--limit N]

RAW_DIR holds gene_list.txt, dataset_list.txt, expression/ and meta/ (the SGD
features file and SPELL's dataset metadata dump). OUT_DIR receives:

* compendium.npz - expression matrices, parsed correctly (no tab collapsing,
  duplicate dataset entries dropped) and normalized by spell_next.qc;
* spell.sqlite  - gene and dataset metadata, QC records and precomputed single-gene results;
* qc_report.tsv - per-dataset QC summary for curators.
"""
import argparse
import gzip
import os
import re
import sqlite3
import time
from concurrent.futures import ProcessPoolExecutor

from spell_next import go, qc
from spell_next.compendium import Compendium
from spell_next.search import search

TOP_GENES = 100
TOP_DATASETS = 100

SCHEMA = """
CREATE TABLE genes (
    id INTEGER PRIMARY KEY,         -- SPELL gene ID
    systematic TEXT NOT NULL,
    standard TEXT,
    sgdid TEXT,
    qualifier TEXT,
    description TEXT
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
    has_metadata INTEGER NOT NULL
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


def read_sgd_features(path):
    """Return {systematic name: (sgdid, standard, aliases, qualifier, description)} for ORFs etc."""
    out = {}
    with open(path, errors='replace') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 16 or not f[3] or f[1] in ('CDS', 'intron', 'five_prime_UTR_intron'):
                continue
            aliases = [a for a in f[5].split('|') if a]
            out.setdefault(f[3], (f[0], f[4] or None, aliases, f[2] or None, f[15] or None))
    return out


def read_dataset_dump(path):
    """Return {filename: dict} from SPELL's SQLdump*_fordownloadfiles.txt export."""
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
                'tags': f[10],
            }
    return out


def write_metadata(db, comp, sgd, dumps):
    db.executescript(SCHEMA)
    gene_rows, name_rows = [], []
    for gid, sys_name in zip(comp.gene_ids, comp.gene_names):
        sgdid, standard, aliases, qualifier, desc = sgd.get(sys_name, (None, None, [], None, None))
        gene_rows.append((gid, sys_name, standard, sgdid, qualifier, desc))
        name_rows.append((sys_name.upper(), gid, 'systematic'))
        if standard:
            name_rows.append((standard.upper(), gid, 'standard'))
        name_rows.extend((a.upper(), gid, 'alias') for a in aliases)
    db.executemany('INSERT INTO genes VALUES (?,?,?,?,?,?)', gene_rows)
    db.executemany('INSERT INTO gene_names VALUES (?,?,?)', name_rows)

    for ds in comp.datasets:
        meta = dumps.get(ds.filename)
        geo = re.match(r'(GSE\d+)', ds.filename)
        db.execute('INSERT INTO datasets VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)', (
            ds.id, ds.filename, geo.group(1) if geo else None,
            meta and meta['pubmed_id'],
            meta['title'] if meta else (geo.group(1) if geo else ds.filename),
            meta and meta['description'], meta and meta['first_author'], meta and meta['authors'],
            meta and meta['journal'], meta and meta['year'], meta and meta['tags'],
            ds.num_conds, len(ds.gene_index), int(meta is not None)))


def _qc_one(values, centre=False):
    diag = qc.diagnose(values)
    scale = 'empty' if diag['empty'] else qc.classify_scale(diag)
    z_before = qc.pair_stats(values)[0]
    new_values, actions, warnings = qc.normalize(values, centre=centre)
    # Report condition shifts on the corrected scale (log units), before any centring.
    rescaled = new_values if not centre else qc.normalize(values)[0]
    spread = qc.diagnose(rescaled).get('cond_median_spread')
    z_after, z_std = qc.pair_stats(new_values)
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


def write_go_tables(db, meta_dir):
    """Store GO terms and propagated SGD annotations for the compendium's genes."""
    obo = os.path.join(meta_dir, 'go-basic.obo')
    gaf = os.path.join(meta_dir, 'gene_association.sgd.gaf.gz')
    if not (os.path.exists(obo) and os.path.exists(gaf)):
        print('go: go-basic.obo or gene_association.sgd.gaf.gz missing, skipped')
        return
    terms = go.read_obo(obo)
    sgdid_to_gene = {s: g for g, s in db.execute('SELECT id, sgdid FROM genes WHERE sgdid IS NOT NULL')}
    gene_terms = go.propagate(go.read_gaf(gaf, sgdid_to_gene, terms), go.ancestors(terms))
    used = set().union(*gene_terms.values()) if gene_terms else set()
    db.executemany('INSERT INTO go_terms VALUES (?,?,?)', [(t, terms[t][0], terms[t][1]) for t in sorted(used)])
    db.executemany('INSERT INTO gene_go VALUES (?,?)', [(g, t) for g, ts in gene_terms.items() for t in ts])
    with open(obo) as fh:
        version = next((ln.split(': ', 1)[1].strip() for ln in fh if ln.startswith('data-version')), 'unknown')
    with gzip.open(gaf, 'rt') as fh:
        gaf_date = next((ln.split(':', 1)[1].strip() for ln in fh if ln.startswith('!date-generated')
                         or ln.startswith('!Generated')), 'unknown')
    db.executemany('INSERT OR REPLACE INTO build_info VALUES (?,?)',
                   [('go_ontology', version), ('go_annotations', gaf_date)])
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


def precompute(db, npz_path, gene_ids, workers):
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('raw_dir')
    ap.add_argument('out_dir')
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--limit', type=int, help='precompute only the first N genes (for testing)')
    ap.add_argument('--no-qc', action='store_true', help='skip normalization (use values as published)')
    ap.add_argument('--centre', action='store_true',
                    help='also centre each condition on its median (off by default: lowers the GO benchmark)')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    npz_path = os.path.join(args.out_dir, 'compendium.npz')
    db_path = os.path.join(args.out_dir, 'spell.sqlite')

    t0 = time.time()
    comp = Compendium.from_spell_files(args.raw_dir, faithful=False, dedupe=True)
    print('compendium: %d datasets, %d genes (%.0fs)' % (len(comp.datasets), len(comp.gene_ids), time.time() - t0))

    if os.path.exists(db_path):
        os.remove(db_path)
    db = sqlite3.connect(db_path)
    meta_dir = os.path.join(args.raw_dir, 'meta')
    sgd = read_sgd_features(os.path.join(meta_dir, 'SGD_features.tab'))
    dump = next(f for f in sorted(os.listdir(meta_dir)) if f.startswith('SQLdump'))
    dumps = read_dataset_dump(os.path.join(meta_dir, dump))
    write_metadata(db, comp, sgd, dumps)
    db.commit()
    n_meta = db.execute('SELECT SUM(has_metadata), COUNT(*) FROM datasets').fetchone()
    print('metadata: %d/%d datasets have publication metadata' % n_meta)

    write_go_tables(db, meta_dir)

    if not args.no_qc:
        apply_qc(db, comp, args.workers, centre=args.centre)
        write_qc_report(db, os.path.join(args.out_dir, 'qc_report.tsv'))
    comp.save(npz_path)

    gene_ids = comp.gene_ids[:args.limit] if args.limit else comp.gene_ids
    print('precomputing %d genes with %d workers...' % (len(gene_ids), args.workers))
    precompute(db, npz_path, gene_ids, args.workers)
    db.execute('VACUUM')
    db.close()
    print('done in %.0fs: %s' % (time.time() - t0, db_path))


if __name__ == '__main__':
    main()
