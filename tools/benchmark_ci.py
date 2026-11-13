"""Paired bootstrap CIs for the GO benchmark, relative to the first build.

Usage: benchmark_ci.py GAF_GZ BASE_BUILD OTHER_BUILD ... [--faithful COMPENDIUM_NPZ]
"""
import sqlite3
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import benchmark_go as B


def main(args):
    faithful = None
    if '--faithful' in args:
        i = args.index('--faithful')
        faithful = args[i + 1]
        args = args[:i] + args[i + 2:]
    gaf, builds = args[0], args[1:]
    db = sqlite3.connect(builds[0] + '/spell.sqlite')
    terms = B.go_terms(gaf, {s: g for g, s in db.execute('SELECT id, sgdid FROM genes WHERE sgdid IS NOT NULL')})
    genes = sorted(terms)

    def per_gene(hits):
        return np.array([sum(bool(terms[g] & terms.get(h, set())) for h in hits[g]) / len(hits[g]) for g in genes])

    res = {b: per_gene(B.hits_from_db(b + '/spell.sqlite')) for b in builds}
    if faithful:
        for prec in ('float32', 'float64'):
            with ProcessPoolExecutor(16, initializer=B._init, initargs=(faithful, prec)) as pool:
                res['old parsing, %s' % prec] = per_gene(dict(pool.map(B._faithful, genes, chunksize=16)))

    rng = np.random.default_rng(1)
    idx = rng.integers(0, len(genes), size=(2000, len(genes)))
    base = res[builds[0]]
    print('difference in GO benchmark vs %s (%d genes, 2000 paired bootstrap resamples)' % (builds[0], len(genes)))
    for name, v in res.items():
        if name == builds[0]:
            continue
        d = v - base
        lo, hi = np.quantile(d[idx].mean(axis=1), [.025, .975])
        print('  %-24s %+.4f  (95%% CI %+.4f to %+.4f)%s' % (
            name, d.mean(), lo, hi, '' if lo < 0 < hi else '  significant'))


if __name__ == '__main__':
    main(sys.argv[1:])
