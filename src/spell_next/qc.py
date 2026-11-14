"""Quality control and normalization of expression datasets.

Problems seen in the SPELL compendium:

* some datasets hold unlogged values (raw intensities or counts), so a few highly
  expressed genes dominate every correlation;
* some datasets mix arrays or batches on different scales, so whole conditions are
  shifted up or down and every gene looks correlated with every other.

``normalize`` log-transforms clearly unlogged datasets and reports condition shifts for
curators. Centring each condition on its median is available (``centre=True``) but off
by default: on a genome-wide GO benchmark it made results significantly worse, because a
condition-wide shift is often real biology (e.g. a strong stress response), not an
artifact. ``pair_stats`` recomputes a dataset's z_mean/z_std (mean and standard
deviation of Fisher z over all gene pairs), which SPELL's scores are standardized by.
"""
import numpy as np

from spell_next.search import centre_rows

UNLOGGED_MAX = 100.0       # log-scale expression values never get this large
SHIFT_WARN = 1.0           # spread of condition medians (log units) worth reporting
FISHER_LIMIT = 0.999999


def pairwise_pearson(values, rows=None):
    """Pearson correlation of the given rows (default: all) with every row, using
    pairwise-complete conditions (NaN = missing). Returns a len(rows) x n matrix."""
    mask = ~np.isnan(values)
    m = mask.astype(np.float64)
    x = np.where(mask, centre_rows(values), 0.0)     # centred: stable for near-constant genes
    rows = slice(None) if rows is None else rows
    mi, xi = m[rows], x[rows]
    n = mi @ m.T
    sx = xi @ m.T                     # sum of row i's values where row j is also present
    sy = mi @ x.T                     # sum of row j's values where row i is also present
    sxx = (xi * xi) @ m.T
    syy = mi @ (x * x).T
    sxy = xi @ x.T
    with np.errstate(divide='ignore', invalid='ignore'):
        r = (sxy - sx * sy / n) / np.sqrt((sxx - sx * sx / n) * (syy - sy * sy / n))
    r[n < 3] = np.nan
    return np.clip(r, -1.0, 1.0)


def pair_stats(values, max_genes=None, seed=0, block=1024):
    """Mean and population SD of Fisher z over all gene pairs (as Spell++ writes .sinfo).

    Rows are processed in blocks so memory stays near block x genes rather than genes x genes.
    """
    if max_genes is not None and values.shape[0] > max_genes:
        rows = np.random.default_rng(seed).choice(values.shape[0], max_genes, replace=False)
        values = values[rows]
    n_genes = values.shape[0]
    total = total_sq = 0.0
    count = 0
    for i0 in range(0, n_genes, block):
        i1 = min(n_genes, i0 + block)
        r = pairwise_pearson(values, slice(i0, i1))
        upper = np.arange(n_genes)[None, :] > np.arange(i0, i1)[:, None]      # pairs j > i only
        r = np.clip(r[upper], -FISHER_LIMIT, FISHER_LIMIT)
        z = 0.5 * np.log((1 + r) / (1 - r))
        z = z[~np.isnan(z)]
        total += float(z.sum())
        total_sq += float((z * z).sum())
        count += z.size
    if count == 0:
        return float('nan'), float('nan')
    mean = total / count
    return mean, float(np.sqrt(max(0.0, total_sq / count - mean * mean)))


def diagnose(values):
    """Summary statistics used to decide how to normalize a dataset."""
    present = values[~np.isnan(values)]
    if present.size == 0:
        return {'empty': True}
    with np.errstate(all='ignore'):
        cond_medians = np.nanmedian(values, axis=0)
    cond_medians = cond_medians[~np.isnan(cond_medians)]
    return {
        'empty': False,
        'min': float(present.min()),
        'median': float(np.median(present)),
        'max': float(present.max()),
        'frac_missing': float(np.isnan(values).mean()),
        'frac_negative': float((present < 0).mean()),
        'q01': float(np.quantile(present, 0.01)),
        'cond_median_spread': float(cond_medians.max() - cond_medians.min()) if cond_medians.size else 0.0,
    }


def classify_scale(d):
    """Return 'log', 'ratio' (unlogged ratios), 'intensity' (unlogged values) or 'unclear'."""
    if d['frac_negative'] > 0.01:
        # Negative values mean log data, possibly offset from 0; an offset does not change correlations.
        return 'unclear' if d['max'] > UNLOGGED_MAX else 'log'
    if 0.5 < d['median'] < 2 and d['max'] > 4:
        # Unlogged ratios straddle 1; all-positive values that never drop below 1 are something else.
        return 'ratio' if d['q01'] < 0.8 else 'unclear'
    if d['max'] > UNLOGGED_MAX:
        return 'intensity'
    return 'log'


def normalize(values, centre=False):
    """Return (normalized values, list of actions taken, list of warnings).

    ``centre=True`` also subtracts each condition's median.
    """
    actions, warnings = [], []
    d = diagnose(values)
    if d['empty']:
        return values, actions, ['no values']
    out = values.astype(np.float64, copy=True)

    scale = classify_scale(d)
    if scale == 'ratio':
        with np.errstate(divide='ignore', invalid='ignore'):
            out = np.where(out > 0, np.log2(out), np.nan)
        actions.append('log2 transform of unlogged ratios (median %.2f)' % d['median'])
    elif scale == 'intensity':
        out = np.log2(np.clip(out, 0, None) + 1)
        actions.append('log2(x+1) transform of unlogged values (median %.0f, max %.0f)' % (d['median'], d['max']))
    elif scale == 'unclear':
        warnings.append('unclear scale (median %.2f, max %.0f, %.0f%% negative): left as is'
                        % (d['median'], d['max'], 100 * d['frac_negative']))

    with np.errstate(all='ignore'):
        medians = np.nanmedian(out, axis=0)
    spread = float(np.nanmax(medians) - np.nanmin(medians)) if np.any(~np.isnan(medians)) else 0.0
    if centre:
        out -= np.where(np.isnan(medians), 0.0, medians)
        actions.append('condition medians centred (spread before %.2f)' % spread)
    if spread > SHIFT_WARN:
        warnings.append('conditions shifted by up to %.2f log units%s' % (spread, ' before centring' if centre else ''))
    return out, actions, warnings
