"""SPELL query-driven co-expression search, ported from Spell++ (Searcher.cpp).

Production SPELL runs with the ZSCORE2 distance and no orthology, so that is what is
implemented here:

* score(x, y) in a dataset = z = (fisherZ(pearson(x, y)) - z_mean) / z_std, using the
  dataset's precomputed z_mean/z_std; z < 1 counts as 0, otherwise z squared.
* Multi-gene query: a dataset's weight is the mean score over pairs of query genes.
  Each gene's score is the weighted sum over datasets of its mean score against the
  query genes, divided by 100.
* Single-gene query: datasets are weighted by the variance of the query gene, but
  that weight only orders the datasets; gene scores use weight 1. The displayed
  dataset order comes from a second pass with the query gene plus its top 2 hits.

``precision='float32'`` reproduces Spell++ arithmetic: single-precision values and
running sums, and the one-pass Pearson formula, which loses accuracy badly for some
gene pairs. ``precision='float64'`` computes the same quantities accurately, on values
centred per gene (the one-pass formula is unstable for near-constant genes otherwise).
"""
import math
from dataclasses import dataclass

import numpy as np

MIN_CONDITIONS = 5


@dataclass
class SearchResult:
    """Ranked datasets and genes, as returned by the Spell++ engine."""

    dataset_ids: list           # ordered by final weight, highest first
    dataset_weights: list       # percent of total weight (None for single-gene queries)
    gene_ids: list              # ordered by score, highest first
    gene_scores: list           # sqrt(score), as SPELL displays it
    raw_weights: list = None    # unnormalized final weight of each dataset in dataset_ids

    def to_engine_lines(self, max_genes=None):
        """Format like Spell++ SearchResult::toString() (two lines)."""
        if self.dataset_weights is None:
            l1 = ','.join('%d:--' % d for d in self.dataset_ids)
        else:
            l1 = ','.join('%d:%.1f%%' % (d, w) for d, w in zip(self.dataset_ids, self.dataset_weights))
        n = len(self.gene_ids) if max_genes is None else max_genes
        l2 = ','.join('%d:%.1f' % (g, s) for g, s in zip(self.gene_ids[:n], self.gene_scores[:n]))
        return l1, l2


class _Arith:
    """Arithmetic in one floating-point type; float32 sums run sequentially like C++ loops."""

    def __init__(self, precision):
        if precision not in ('float32', 'float64'):
            raise ValueError('precision must be float32 or float64')
        self.dtype = np.float32 if precision == 'float32' else np.float64

    def f(self, x):
        return self.dtype(x)

    def rowsum(self, a):
        if self.dtype is np.float64:
            return a.sum(axis=1)
        if a.shape[1] == 0:
            return np.zeros(a.shape[0], dtype=np.float32)
        return np.cumsum(a, axis=1, dtype=np.float32)[:, -1]

    def seqsum(self, v):
        if self.dtype is np.float64:
            return float(np.sum(v))
        return np.cumsum(v, dtype=np.float32)[-1] if len(v) else np.float32(0)

    def values(self, ds):
        if self.dtype is np.float64:
            # Centre each gene on its mean: correlations are unchanged mathematically, but the
            # one-pass formula no longer loses all precision for near-constant genes.
            if getattr(ds, '_centred', None) is None:
                ds._centred = centre_rows(ds.values)
            return ds._centred
        if getattr(ds, '_values32', None) is None:
            ds._values32 = ds.values.astype(np.float32)
        return ds._values32


def centre_rows(values):
    """Subtract each row's mean over its present values (rows with no values stay NaN)."""
    present = ~np.isnan(values)
    n = present.sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = np.where(present, values, 0.0).sum(axis=1) / n
    return values - np.where(n > 0, mean, 0.0)[:, None]


def pearson_with_row(values, row, ar):
    """Pearson correlation of one row against every row, using pairwise-complete values."""
    y = values[row]
    valid = ~np.isnan(values) & ~np.isnan(y)
    n = valid.sum(axis=1).astype(ar.dtype)
    zero = ar.f(0)
    xv = np.where(valid, values, zero)
    yv = np.where(valid, y, zero)
    sx, sy = ar.rowsum(xv), ar.rowsum(yv)
    sxx, syy, sxy = ar.rowsum(xv * xv), ar.rowsum(yv * yv), ar.rowsum(xv * yv)
    with np.errstate(divide='ignore', invalid='ignore'):
        rho = (sxy - (sx * sy) / n) / np.sqrt((sxx - (sx * sx) / n) * (syy - (sy * sy) / n))
    one = ar.f(1)
    rho = np.where(rho > one, one, np.where(rho < -one, -one, rho))
    rho[n < 3] = np.nan
    return rho


def fisher_z(r, ar):
    lim = ar.f(0.999999)
    if ar.dtype is np.float64:
        # Cap |r| at 0.999999 like Spell++ intends. Spell++ only caps |r| >= 1, which in double
        # precision lets values a rounding step below 1 through (z ~ 19 instead of 7.3), so a
        # single near-constant gene pair could outweigh hundreds of datasets.
        r = np.clip(r, -lim, lim)
    else:
        r = np.where(r <= -1, -lim, np.where(r >= 1, lim, r)).astype(ar.dtype)
    with np.errstate(divide='ignore', invalid='ignore'):
        return ar.f(0.5) * np.log((ar.f(1) + r) / (ar.f(1) - r))


def zscore2(ds, rho, ar):
    """ZSCORE2 distance: standardized Fisher z, zeroed below 1 and squared above."""
    z = (fisher_z(rho, ar) - ar.f(ds.z_mean)) / ar.f(ds.z_std)
    return np.where(np.isnan(z), np.nan, np.where(z < 1, ar.f(0), z * z)).astype(ar.dtype)


def gene_variance(ds, gene, ar):
    """Spell++'s single-gene dataset weight (its stddev() actually returns the variance)."""
    row = ds.row_of.get(gene)
    if row is None:
        return math.nan
    x = ar.values(ds)[row]
    x = x[~np.isnan(x)]
    if x.size == 0:
        return math.nan
    mean = ar.seqsum(x) / ar.f(x.size)
    return ar.seqsum((x - mean) * (x - mean)) / ar.f(x.size)


def pair_score(ds, gi, gj, ar):
    ri, rj = ds.row_of.get(gi), ds.row_of.get(gj)
    if ri is None or rj is None:
        return math.nan
    rho = pearson_with_row(ar.values(ds)[[ri, rj]], 0, ar)[1:]
    return zscore2(ds, rho, ar)[0]


def dataset_weights(datasets, query, ar):
    """Return (kept datasets, their weights, num_weighted) as Searcher::querySearch does."""
    weights, num_weighted = [], 0
    for ds in datasets:
        if len(query) == 1:
            w = gene_variance(ds, query[0], ar)
            if not math.isnan(w):
                num_weighted += 1
        else:
            w, pairs = ar.f(0), 0
            for i in range(len(query)):
                for j in range(i + 1, len(query)):
                    s = pair_score(ds, query[i], query[j], ar)
                    if not math.isnan(s):
                        pairs += 1
                        w = ar.f(w + s)
            if pairs > 0 and w > 0:
                w = ar.f(w / ar.f(pairs))
                num_weighted += 1
        weights.append(w)

    kept, kept_w = [], []
    for ds, w in zip(datasets, weights):
        if not math.isnan(w) and ds.num_conds >= MIN_CONDITIONS:
            kept.append(ds)
            kept_w.append(w)
        elif w > 0:
            num_weighted -= 1

    if num_weighted != 0 and num_weighted < 3:
        kept_w = [ar.f(1)] * len(kept)
        num_weighted = len(kept)
    return kept, kept_w, num_weighted


def score_genes(num_genes, datasets, weights, num_weighted, query, ar):
    """Score every gene in the compendium; returns raw scores (before sqrt)."""
    single = len(query) == 1
    total = np.zeros(num_genes, dtype=ar.dtype)
    used = np.zeros(num_genes, dtype=np.int64)
    for ds, w in zip(datasets, weights):
        if not w > 0:
            continue
        values = ar.values(ds)
        acc = np.zeros(len(ds.gene_index), dtype=ar.dtype)
        cnt = np.zeros(len(ds.gene_index), dtype=np.int64)
        for q in query:
            row = ds.row_of.get(q)
            if row is None:
                continue
            vals = zscore2(ds, pearson_with_row(values, row, ar), ar)
            ok = ~np.isnan(vals)
            acc[ok] += vals[ok]
            cnt[ok] += 1
        has = cnt > 0
        wgt = ar.f(1) if single else ar.f(w)
        total[ds.gene_index[has]] += wgt * (acc[has] / cnt[has].astype(ar.dtype))
        used[ds.gene_index[has]] += 1
    total /= ar.f(100)
    threshold = min(float(num_weighted), max(0.1 * len(datasets), 4.0))
    total[used < threshold] = 0
    return total


def search(compendium, query_ids, dataset_ids=None, max_genes=100, precision='float64'):
    """Run a SPELL search. ``query_ids`` are SPELL gene IDs (gene_list.txt)."""
    ar = _Arith(precision)
    query = [compendium.index_of_id[g] for g in query_ids]
    wanted = None if dataset_ids is None else set(dataset_ids)
    datasets = [d for d in compendium.datasets if wanted is None or d.id in wanted]

    kept, weights, num_weighted = dataset_weights(datasets, query, ar)
    scores = score_genes(len(compendium.gene_ids), kept, weights, num_weighted, query, ar)
    order = np.argsort(-scores, kind='stable')

    final_weights = list(weights)
    if len(query) == 1:
        pass2_query = query + [int(order[1]), int(order[2])]
        p2_kept, p2_weights, _ = dataset_weights(datasets, pass2_query, ar)
        p2 = {ds.id: w for ds, w in zip(p2_kept, p2_weights)}
        final_weights = [p2.get(ds.id, w) for ds, w in zip(kept, weights)]

    total_w = ar.seqsum(np.array(final_weights, dtype=ar.dtype))
    ds_order = sorted(range(len(kept)), key=lambda k: -final_weights[k])
    # Spell++ returns the requested number of genes plus the query genes themselves.
    top = order if max_genes is None else order[:max_genes + len(query)]
    return SearchResult(
        dataset_ids=[kept[k].id for k in ds_order],
        dataset_weights=None if len(query) == 1 else
        [float(ar.f(final_weights[k] / total_w * ar.f(100))) for k in ds_order],
        gene_ids=[compendium.gene_ids[i] for i in top],
        gene_scores=[float(np.sqrt(scores[i])) for i in top],
        raw_weights=[float(final_weights[k]) for k in ds_order],
    )
