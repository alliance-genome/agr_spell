"""Benchmark builds on known co-regulated gene modules and on dataset dominance.

Usage: benchmark.py BUILD_DIR [BUILD_DIR ...] [--faithful COMPENDIUM_NPZ]

For each module gene, precision@20 = fraction of its top 20 hits in the same module.
"""
import re
import sqlite3
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

MODULES = {
    'cytosolic ribosome': re.compile(r'^RP[LS]\d+[AB]?$|^RPP[0-2][AB]?$'),
    'proteasome': re.compile(r'^(PRE\d+|PUP[1-3]|SCL1|RPN\d+|RPT[1-6])$'),
    'histones': re.compile(r'^(HHT[12]|HHF[12]|HTA[12]|HTB[12])$'),
}
K = 20


def modules(db):
    out = {}
    for gid, std in db.execute('SELECT id, standard FROM genes WHERE standard IS NOT NULL'):
        for name, pattern in MODULES.items():
            if pattern.match(std):
                out.setdefault(name, set()).add(gid)
    return out


def precision_from_db(db, mods):
    res = {}
    for name, genes in mods.items():
        vals = []
        for g in genes:
            hits = [h for (h,) in db.execute('SELECT hit_id FROM coexpression WHERE gene_id=? AND rank<=? '
                                             'ORDER BY rank', (g, K))]
            if hits:
                vals.append(sum(h in genes for h in hits) / len(hits))
        res[name] = (float(np.mean(vals)), len(vals))
    return res


_comp = None


def _init(npz):
    global _comp
    from spell_next.compendium import Compendium
    _comp = Compendium.load(npz)


def _faithful_hits(gene_id):
    from spell_next.search import search
    r = search(_comp, [gene_id], max_genes=K, precision='float32')
    return gene_id, [g for g in r.gene_ids if g != gene_id][:K], r.dataset_ids[:10]


def main(args):
    faithful = None
    if '--faithful' in args:
        i = args.index('--faithful')
        faithful = args[i + 1]
        args = args[:i] + args[i + 2:]
    builds = args
    db0 = sqlite3.connect(builds[0] + '/spell.sqlite')
    mods = modules(db0)
    print('modules: ' + ', '.join('%s (%d genes)' % (k, len(v)) for k, v in mods.items()))

    rows = {}
    if faithful:
        genes = sorted(set().union(*mods.values()))
        with ProcessPoolExecutor(16, initializer=_init, initargs=(faithful,)) as pool:
            hits = {g: h for g, h, _ in pool.map(_faithful_hits, genes)}
        rows['old SPELL (faithful)'] = {
            name: (float(np.mean([sum(h in gs for h in hits[g]) / len(hits[g]) for g in gs if hits[g]])), len(gs))
            for name, gs in mods.items()}
    for b in builds:
        rows[b] = precision_from_db(sqlite3.connect(b + '/spell.sqlite'), mods)

    print('\nprecision@%d within module (mean over module genes)' % K)
    print('%-28s' % 'build' + ''.join('%22s' % m for m in mods))
    for b, res in rows.items():
        print('%-28s' % b + ''.join('%22.3f' % res[m][0] for m in mods))

    print('\ndataset dominance: how many genes list a dataset among their top 10 datasets')
    for b in builds:
        db = sqlite3.connect(b + '/spell.sqlite')
        counts = [c for (c,) in db.execute('SELECT COUNT(*) FROM coexpression_datasets WHERE rank<=10 '
                                           'GROUP BY dataset_id ORDER BY 1 DESC')]
        top = db.execute('SELECT d.filename, COUNT(*) FROM coexpression_datasets c JOIN datasets d '
                         'ON d.id=c.dataset_id WHERE c.rank<=10 GROUP BY c.dataset_id '
                         'ORDER BY 2 DESC LIMIT 3').fetchall()
        print('  %-24s max %4d, top-5 datasets share %4.1f%% of all top-10 slots; most frequent: %s' % (
            b, counts[0], 100 * sum(counts[:5]) / sum(counts), ', '.join('%s (%d)' % t for t in top)))


if __name__ == '__main__':
    main(sys.argv[1:])
