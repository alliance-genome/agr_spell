"""Genome-wide benchmark: do a gene's top hits share specific GO biological-process terms?

Usage: benchmark_go.py GAF_GZ BUILD_DIR [BUILD_DIR ...] [--faithful COMPENDIUM_NPZ]

A term is "specific" if 5-100 genes are directly annotated to it (NOT annotations skipped).
For each annotated gene: share of its top 20 hits annotated to at least one of its terms.
Also reports how much rankings differ between the first two builds.
"""
import gzip
import sqlite3
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np

K = 20


def go_terms(gaf, sgdid_to_gene):
    by_term = defaultdict(set)
    with gzip.open(gaf, 'rt') as fh:
        for line in fh:
            if line.startswith('!'):
                continue
            f = line.split('\t')
            if len(f) < 9 or f[8] != 'P' or 'NOT' in f[3]:
                continue
            gene = sgdid_to_gene.get(f[1])
            if gene is not None:
                by_term[f[4]].add(gene)
    terms = defaultdict(set)
    for term, genes in by_term.items():
        if 5 <= len(genes) <= 100:
            for g in genes:
                terms[g].add(term)
    return terms


def score(hits_of, terms):
    vals = []
    for g, hits in hits_of.items():
        if g in terms and hits:
            vals.append(sum(bool(terms[g] & terms.get(h, set())) for h in hits) / len(hits))
    return float(np.mean(vals)), len(vals)


def hits_from_db(path):
    out = defaultdict(list)
    for g, h in sqlite3.connect(path).execute(
            'SELECT gene_id, hit_id FROM coexpression WHERE rank <= ? ORDER BY gene_id, rank', (K,)):
        out[g].append(h)
    return out


_comp = None


_precision = 'float32'


def _init(npz, precision='float32'):
    global _comp, _precision
    _precision = precision
    from spell_next.compendium import Compendium
    _comp = Compendium.load(npz)


def _faithful(gene_id):
    from spell_next.search import search
    r = search(_comp, [gene_id], max_genes=K, precision=_precision)
    return gene_id, [g for g in r.gene_ids if g != gene_id][:K]


def main(args):
    faithful = None
    if '--faithful' in args:
        i = args.index('--faithful')
        faithful = args[i + 1]
        args = args[:i] + args[i + 2:]
    gaf, builds = args[0], args[1:]
    db = sqlite3.connect(builds[0] + '/spell.sqlite')
    sgdid_to_gene = {s: g for g, s in db.execute('SELECT id, sgdid FROM genes WHERE sgdid IS NOT NULL')}
    terms = go_terms(gaf, sgdid_to_gene)
    print('genes with a specific GO BP term: %d' % len(terms))

    results = {}
    if faithful:
        genes = [g for (g,) in db.execute('SELECT id FROM genes')]
        with ProcessPoolExecutor(16, initializer=_init, initargs=(faithful,)) as pool:
            results['old SPELL (faithful)'] = dict(pool.map(_faithful, genes, chunksize=16))
    for b in builds:
        results[b] = hits_from_db(b + '/spell.sqlite')

    print('\nshare of top-%d hits sharing a specific GO process term (higher is better)' % K)
    for name, hits in results.items():
        mean, n = score(hits, terms)
        print('  %-24s %.4f  (%d genes)' % (name, mean, n))
    rng = np.random.default_rng(0)
    shuffled = {g: list(rng.choice(list(results[builds[0]]), K, replace=False)) for g in results[builds[0]]}
    print('  %-24s %.4f' % ('random hits', score(shuffled, terms)[0]))

    names = list(results)
    for a, b in zip(names, names[1:]):
        ha, hb = results[a], results[b]
        overlap = [len(set(ha[g]) & set(hb[g])) for g in ha if g in hb and ha[g]]
        print('\n%s vs %s: mean top-%d overlap %.1f; genes with <10 shared: %d of %d'
              % (a, b, K, np.mean(overlap), sum(o < 10 for o in overlap), len(overlap)))


if __name__ == '__main__':
    main(sys.argv[1:])
