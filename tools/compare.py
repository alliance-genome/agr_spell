"""Compare spell-next search results with the live Spell++ engine baseline.

Usage: compare.py COMPENDIUM_NPZ QUERY_IDS_TSV BASELINE_TSV [--precision=float32|float64] [--verbose]
"""
import sys
import time
from collections import defaultdict

from spell_next.compendium import Compendium
from spell_next.search import search


def parse_genes(line):
    out = []
    for item in line.split(','):
        if item:
            gid, score = item.rsplit(':', 1)
            out.append((int(gid), score))
    return out


def parse_datasets(line):
    return [tuple(item.split(':', 1)) for item in line.split(',') if item]


def same_up_to_ties(expected, actual):
    """True if both lists hold the same items under each displayed value, in value order."""
    if [v for _, v in expected] != [v for _, v in actual]:
        return False
    groups_e, groups_a = defaultdict(set), defaultdict(set)
    for k, v in expected:
        groups_e[v].add(k)
    for k, v in actual:
        groups_a[v].add(k)
    return groups_e == groups_a


def drop_last_group(items):
    """Drop the trailing tie group, which a top-N cut may split arbitrarily."""
    if not items:
        return items
    last = items[-1][1]
    return [x for x in items if x[1] != last]


def order_consistent(expected_ids, weight_of):
    """True if expected_ids is a valid descending order under our weights (ties free)."""
    try:
        w = [weight_of[d] for d in expected_ids]
    except KeyError:
        return False
    return all(a >= b for a, b in zip(w, w[1:]))


def main(npz, query_file, baseline_file, verbose=False, precision='float32'):
    t0 = time.time()
    comp = Compendium.load(npz)
    print('loaded %d datasets, %d genes in %.1fs' % (len(comp.datasets), len(comp.gene_ids), time.time() - t0))

    queries = {}
    with open(query_file) as fh:
        for line in fh:
            qid, genes = line.rstrip('\n').split('\t')
            queries[qid] = [int(g) for g in genes.split(',')]

    stats = defaultdict(int)
    times = []
    with open(baseline_file) as fh:
        for line in fh:
            qid, l1, l2 = line.rstrip('\n').split('\t')
            kind = 'single' if qid.startswith('s') else 'multi'
            t = time.time()
            result = search(comp, queries[qid], precision=precision)
            times.append(time.time() - t)
            a1, a2 = result.to_engine_lines()

            exp_g, act_g = parse_genes(l2), parse_genes(a2)
            exp_nz = drop_last_group([x for x in exp_g if x[1] != '0.0'])
            act_nz = act_g[:len(exp_nz)]
            gene_ok = same_up_to_ties(exp_nz, act_nz)
            top10 = len({g for g, _ in exp_g[:10]} & {g for g, _ in act_g[:10]})

            exp_d, act_d = parse_datasets(l1), parse_datasets(a1)
            if kind == 'multi':
                ds_ok = same_up_to_ties(exp_d, act_d)
            else:
                weight_of = {str(d): w for d, w in zip(result.dataset_ids, result.raw_weights)}
                ds_ok = order_consistent([d for d, _ in exp_d], weight_of)
            ds_set_ok = {d for d, _ in exp_d} == {d for d, _ in act_d}

            stats[kind + ' queries'] += 1
            stats[kind + ' genes identical'] += gene_ok
            stats[kind + ' top10 identical'] += top10 == 10
            stats[kind + ' datasets identical'] += ds_ok
            stats[kind + ' dataset set identical'] += ds_set_ok
            if verbose and not (gene_ok and ds_ok):
                print('%s genes_ok=%s top10=%d/10 ds_ok=%s ds_set_ok=%s (n_ds exp %d act %d)' % (
                    qid, gene_ok, top10, ds_ok, ds_set_ok, len(exp_d), len(act_d)))
                if not gene_ok:
                    for i, (e, a) in enumerate(zip(exp_nz, act_nz)):
                        if e != a:
                            print('   first gene diff at rank %d: expected %s got %s' % (i, e, a))
                            break
                if not ds_ok:
                    for i, (e, a) in enumerate(zip(exp_d, act_d)):
                        if e != a:
                            print('   first dataset diff at rank %d: expected %s got %s' % (i, e, a))
                            break

    for k in sorted(stats):
        print('%-32s %d' % (k, stats[k]))
    print('search time: median %.2fs, max %.2fs' % (sorted(times)[len(times) // 2], max(times)))


if __name__ == '__main__':
    prec = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--precision=')), 'float32')
    print('precision:', prec)
    main(*sys.argv[1:4], verbose='--verbose' in sys.argv, precision=prec)
