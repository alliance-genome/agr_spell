"""Precompute single-gene SPELL results for every gene at once.

Equivalent to calling ``search(compendium, [gene])`` for each gene, but organised by
dataset: each dataset's gene-by-gene correlations are computed once (in row blocks) and
added to a shared genes x genes score matrix, instead of being recomputed for every
query gene. The second pass (dataset weights from the query gene plus its top two hits)
is vectorised over all genes per dataset. Genes with fewer than 3 usable datasets, where
Spell++ switches to equal weights, are delegated to ``search`` itself.
"""
import math
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from spell_next.compendium import Compendium
from spell_next.qc import pairwise_pearson
from spell_next.search import MIN_CONDITIONS, centre_rows, search

FISHER_LIMIT = 0.999999
BLOCK = 512


def zscore2(r, z_mean, z_std):
    """ZSCORE2 of correlations (as spell_next.search.zscore2, float64); NaN stays NaN."""
    r = np.clip(r, -FISHER_LIMIT, FISHER_LIMIT)
    with np.errstate(invalid='ignore', divide='ignore'):
        z = (0.5 * np.log((1 + r) / (1 - r)) - z_mean) / z_std
        return np.where(np.isnan(z), np.nan, np.where(z < 1, 0.0, z * z))


def gene_variances(values):
    """Population variance of each row over its present values (NaN if none)."""
    present = ~np.isnan(values)
    n = present.sum(axis=1)
    x = np.where(present, values, 0.0)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = x.sum(axis=1) / n
        dev = np.where(present, values - mean[:, None], 0.0)
        var = (dev * dev).sum(axis=1) / n
    var[n == 0] = np.nan
    return var


def pair_pearson(a, b):
    """Row-wise Pearson correlation of a[i] with b[i], using pairwise-complete conditions."""
    valid = ~np.isnan(a) & ~np.isnan(b)
    n = valid.sum(axis=1)
    x, y = np.where(valid, a, 0.0), np.where(valid, b, 0.0)
    sx, sy = x.sum(axis=1), y.sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        r = ((x * y).sum(axis=1) - sx * sy / n) / np.sqrt(
            ((x * x).sum(axis=1) - sx * sx / n) * ((y * y).sum(axis=1) - sy * sy / n))
    r = np.clip(r, -1.0, 1.0)
    r[n < 3] = np.nan
    return r


_comp = None


def _init(npz_path):
    global _comp
    _comp = Compendium.load(npz_path)


def _pass1(ds_indices):
    """Sum of ZSCORE2 scores and count of contributing datasets, for every (query, gene) pair."""
    n = len(_comp.gene_ids)
    total = np.zeros((n, n))
    used = np.zeros((n, n), dtype=np.int32)
    for k in ds_indices:
        ds = _comp.datasets[k]
        if ds.num_conds < MIN_CONDITIONS:
            continue
        gi = np.asarray(ds.gene_index)
        rows = np.nonzero(gene_variances(ds.values) > 0)[0]     # single-gene weight > 0
        for b in range(0, len(rows), BLOCK):
            rb = rows[b:b + BLOCK]
            z2 = zscore2(pairwise_pearson(ds.values, rb), ds.z_mean, ds.z_std)
            ok = ~np.isnan(z2)
            idx = np.ix_(gi[rb], gi)
            total[idx] += np.where(ok, z2, 0.0)
            used[idx] += ok
    return total, used


def _pass2(args):
    """Per dataset and query gene: summed ZSCORE2 over the query's three gene pairs, and pair count."""
    ds_indices, triples = args
    n = len(triples)
    w_sum = np.zeros((len(ds_indices), n))
    pairs = np.zeros((len(ds_indices), n), dtype=np.int8)
    for j, k in enumerate(ds_indices):
        ds = _comp.datasets[k]
        values = centre_rows(ds.values)
        row_of = np.full(len(_comp.gene_ids), -1)
        row_of[np.asarray(ds.gene_index)] = np.arange(len(ds.gene_index))
        rows = row_of[triples]                                   # n x 3, -1 = not measured
        for a, b in ((0, 1), (0, 2), (1, 2)):
            both = (rows[:, a] >= 0) & (rows[:, b] >= 0)
            if not both.any():
                continue
            s = zscore2(pair_pearson(values[rows[both, a]], values[rows[both, b]]), ds.z_mean, ds.z_std)
            ok = ~np.isnan(s)
            idx = np.nonzero(both)[0][ok]
            w_sum[j, idx] += s[ok]
            pairs[j, idx] += 1
    return ds_indices, w_sum, pairs


def _chunks(items, k):
    k = max(1, min(k, len(items)))
    return [items[i::k] for i in range(k)]


def precompute_all(comp, npz_path, top_genes=100, top_datasets=100, workers=None, log=print):
    """Return {gene_id: (hits, datasets)} as build._search_one would for every gene.

    hits = [(gene_id, sqrt score)] (query excluded), datasets = [(dataset_id, percent weight)].
    """
    workers = workers or 1
    n = len(comp.gene_ids)
    datasets = comp.datasets
    n_ds = len(datasets)
    small = np.array([ds.num_conds < MIN_CONDITIONS for ds in datasets])

    # Single-gene dataset weights (variance) per dataset and gene, as Spell++ computes them.
    var = np.full((n_ds, n), np.nan)
    for k, ds in enumerate(datasets):
        var[k, np.asarray(ds.gene_index)] = gene_variances(ds.values)
    present = ~np.isnan(var)
    kept1 = present & ~small[:, None]                                  # datasets searched per gene
    num_weighted = present.sum(axis=0) - ((var > 0) & small[:, None]).sum(axis=0)
    n_kept = kept1.sum(axis=0)
    threshold = np.minimum(num_weighted, np.maximum(0.1 * n_kept, 4.0))
    special = num_weighted < 3                                          # Spell++'s equal-weights rule

    with ProcessPoolExecutor(workers, initializer=_init, initargs=(npz_path,)) as pool:
        total = np.zeros((n, n))
        used = np.zeros((n, n), dtype=np.int32)
        for t, u in pool.map(_pass1, _chunks(list(range(n_ds)), workers)):
            total += t
            used += u
            del t, u
        log('  pass 1 (scores) done')

        scores = total / 100
        del total
        scores[used < threshold[:, None]] = 0.0
        del used
        order = np.argsort(-scores, axis=1, kind='stable')
        triples = np.column_stack([np.arange(n), order[:, 1], order[:, 2]])

        w_sum = np.zeros((n_ds, n))
        pairs = np.zeros((n_ds, n), dtype=np.int8)
        jobs = [(c, triples) for c in _chunks(list(range(n_ds)), workers)]
        for idx, ws, ps in pool.map(_pass2, jobs):
            w_sum[idx] = ws
            pairs[idx] = ps
        log('  pass 2 (dataset weights) done')

    positive = (pairs > 0) & (w_sum > 0)
    with np.errstate(invalid='ignore', divide='ignore'):
        w2 = np.where(positive, w_sum / np.maximum(pairs, 1), w_sum)
    num_weighted2 = (positive & ~small[:, None]).sum(axis=0)
    equal2 = (num_weighted2 != 0) & (num_weighted2 < 3)

    ds_ids = np.array([ds.id for ds in datasets])
    results = {}
    for q in range(n):
        gid = comp.gene_ids[q]
        if special[q]:
            continue
        top = order[q, :top_genes + 1]
        hits = [(comp.gene_ids[g], math.sqrt(scores[q, g])) for g in top if g != q][:top_genes]
        ks = np.nonzero(kept1[:, q])[0]
        w = np.ones(len(ks)) if equal2[q] else w2[ks, q]
        rank = np.argsort(-w, kind='stable')
        tot = w.sum() or 1.0
        dsets = [(int(ds_ids[ks[r]]), float(w[r] / tot * 100)) for r in rank[:top_datasets]]
        results[gid] = (hits, dsets)

    for q in np.nonzero(special)[0]:
        gid = comp.gene_ids[q]
        r = search(comp, [gid], max_genes=top_genes, precision='float64')
        hits = [(g, s) for g, s in zip(r.gene_ids, r.gene_scores) if g != gid][:top_genes]
        tot = sum(r.raw_weights) or 1.0
        results[gid] = (hits, [(d, w / tot * 100) for d, w in zip(r.dataset_ids, r.raw_weights)][:top_datasets])
    log('  %d genes (%d via the per-gene search)' % (n, int(special.sum())))
    return results
