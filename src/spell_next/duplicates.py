"""Find datasets that duplicate each other: the same expression data, the same conditions, or the same GEO series.

Each expression file gets two fingerprints:

* ``data``: a hash of its values per gene (rounded to 3 decimals, rows in any order), so the
  same matrix saved under another name, or with its rows reordered, matches;
* ``conds``: a hash of its condition names (case and spacing ignored), which catches the same
  experiment reprocessed differently (e.g. ratios vs log ratios). Only for files with at least
  5 conditions, since short generic headers ("0 min, 10 min, ...") repeat across studies.

Fingerprints are cached in ``raw/.dataset_fingerprints.json`` by file size and modification
time, so only new or changed files are read.
"""
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from spell_next.compendium import read_pcl

CACHE = '.dataset_fingerprints.json'
MIN_CONDS = 5


def fingerprint(path, name_to_index):
    """{'data': hex, 'conds': hex or None, 'genes': n} for one PCL file (None if it cannot be read)."""
    try:
        conditions, gene_index, values = read_pcl(path, name_to_index, faithful=False)
    except (OSError, ValueError, IndexError):
        return None
    if len(gene_index) == 0:
        return None
    rounded = np.round(values, 3)
    rows = sorted('%d:%s' % (g, ','.join('nan' if np.isnan(v) else '%.3f' % v for v in row))
                  for g, row in zip(np.asarray(gene_index).tolist(), rounded))
    data = hashlib.sha256('\n'.join(rows).encode()).hexdigest()
    names = [' '.join(c.lower().split()) for c in conditions]
    conds = hashlib.sha256('\t'.join(names).encode()).hexdigest() if len(names) >= MIN_CONDS else None
    return {'data': data, 'conds': conds, 'genes': len(gene_index)}


def fingerprints(raw_dir, paths, name_to_index, workers=None):
    """{path: fingerprint} for the given files, reusing raw/.dataset_fingerprints.json where still valid."""
    cache_path = os.path.join(raw_dir, CACHE)
    try:
        with open(cache_path) as fh:
            cache = json.load(fh)
    except (OSError, ValueError):
        cache = {}
    out, todo = {}, []
    for p in paths:
        st = os.stat(p)
        key = os.path.relpath(p, raw_dir)
        hit = cache.get(key)
        if hit and hit['size'] == st.st_size and hit['mtime_ns'] == st.st_mtime_ns:
            out[p] = hit['fp']
        else:
            todo.append(p)
    if todo:
        with ProcessPoolExecutor(workers) as pool:
            for p, fp in zip(todo, pool.map(fingerprint, todo, [name_to_index] * len(todo), chunksize=4)):
                out[p] = fp
                st = os.stat(p)
                cache[os.path.relpath(p, raw_dir)] = {'size': st.st_size, 'mtime_ns': st.st_mtime_ns, 'fp': fp}
        try:
            tmp = cache_path + '.tmp'
            with open(tmp, 'w') as fh:
                json.dump(cache, fh)
            os.replace(tmp, cache_path)
        except OSError:
            pass                                    # the cache only saves time
    return out


def compare(new, existing, geo_of):
    """Matches of each new dataset against existing ones and earlier new ones.

    new, existing: lists of (label, path, fingerprint); geo_of(label) gives its GEO series or None.
    Returns {label: [(kind, other_label)]} with kind 'same data', 'same conditions' or 'same GEO series'.
    """
    found = {}
    seen = list(existing)
    for label, _, fp in new:
        hits = []
        for other, _, ofp in seen:
            if fp and ofp and fp['data'] == ofp['data']:
                hits.append(('same data', other))
            elif fp and ofp and fp['conds'] and fp['conds'] == ofp['conds']:
                hits.append(('same conditions', other))
            elif geo_of(label) and geo_of(label) == geo_of(other):
                hits.append(('same GEO series', other))
        if hits:
            found[label] = hits
        seen.append((label, None, fp))
    return found


def groups(items, geo_of):
    """Groups of existing datasets that duplicate each other: [(kind, [labels])], most serious first."""
    out = []
    for kind, key in (('same data', lambda fp, label: fp and fp['data']),
                      ('same conditions', lambda fp, label: fp and fp['conds']),
                      ('same GEO series', lambda fp, label: geo_of(label))):
        by = {}
        for label, _, fp in items:
            k = key(fp, label)
            if k:
                by.setdefault(k, []).append(label)
        out += [(kind, labels) for labels in by.values() if len(labels) > 1]
    return out
