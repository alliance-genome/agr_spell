"""Check a MOD's raw inputs before building: config, lists, PCL files, metadata, QC preview.

Errors stop a rebuild (missing or unreadable files, duplicate IDs, datasets with no
recognised genes); warnings are reported but do not (few conditions, low gene-name
match, odd value scale, missing metadata, unknown assay, RNA-seq values that are not
log-scale).
"""
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

import numpy as np

from spell_next import modconfig, qc
from spell_next.compendium import read_dataset_list, read_gene_list, read_pcl, tokenize

MIN_CONDITIONS = 5
MIN_MATCH = 0.5            # warn when fewer than half of a file's rows match known genes
MAX_MISSING = 0.5
COUNT_LIKE = 0.95          # RNA-seq values this often whole numbers (and unlogged) look like raw read counts


@dataclass
class Issue:
    level: str             # error | warning | info
    where: str
    message: str


@dataclass
class Report:
    issues: list = field(default_factory=list)
    datasets: list = field(default_factory=list)     # one dict per dataset
    next_dataset_id: int = 1

    def add(self, level, where, message):
        self.issues.append(Issue(level, where, message))

    @property
    def errors(self):
        return [i for i in self.issues if i.level == 'error']

    @property
    def warnings(self):
        return [i for i in self.issues if i.level == 'warning']

    def write_tsv(self, path):
        cols = ['dataset_id', 'filename', 'conditions', 'rows', 'genes_matched', 'match_rate', 'scale',
                'cond_median_spread', 'frac_missing', 'issues']
        with open(path, 'w') as out:
            out.write('\t'.join(cols) + '\n')
            for d in self.datasets:
                out.write('\t'.join('' if d.get(c) is None else
                                    ('%.3f' % d[c] if isinstance(d[c], float) else str(d[c])) for c in cols) + '\n')


def check_pcl(path, name_to_index):
    """Parse one PCL file the way the build does and describe it (run in a worker process)."""
    out = {'issues': []}
    try:
        with open(path, errors='replace') as fh:
            lines = [ln.rstrip('\r\n') for ln in fh]
        header = tokenize(lines[0], faithful=False) if lines else []
        offset = 3 if len(header) > 2 and header[2].upper() == 'GWEIGHT' else 2
        conditions = header[offset:]
        data_rows = [ln for ln in lines[1:] if ln and not ln.upper().startswith('EWEIGHT')]
        _, gene_index, values = read_pcl(path, name_to_index, faithful=False)
    except (OSError, ValueError, IndexError) as err:
        out['issues'].append(('error', 'cannot be read: %s' % err))
        return out
    out.update(conditions=len(conditions), rows=len(data_rows), genes_matched=len(gene_index),
               match_rate=len(gene_index) / len(data_rows) if data_rows else 0.0)
    if len(gene_index) == 0:
        out['issues'].append(('error', 'no rows match a gene in gene_list.txt'))
        return out
    if out['match_rate'] < MIN_MATCH:
        out['issues'].append(('warning', 'only %.0f%% of rows match known genes' % (100 * out['match_rate'])))
    if len(conditions) < MIN_CONDITIONS:
        out['issues'].append(('warning', '%d conditions: datasets with fewer than %d are not searched'
                              % (len(conditions), MIN_CONDITIONS)))
    if any(not c.strip() for c in conditions):
        out['issues'].append(('warning', 'empty condition name(s)'))
    if len(set(conditions)) < len(conditions):
        out['issues'].append(('warning', 'duplicate condition names'))
    diag = qc.diagnose(values)
    if diag['empty']:
        out['issues'].append(('error', 'no numeric values'))
        return out
    scale = qc.classify_scale(diag)
    rescaled, actions, warnings = qc.normalize(values)
    present = values[~np.isnan(values)]
    out.update(scale=scale, frac_missing=diag['frac_missing'],
               cond_median_spread=qc.diagnose(rescaled).get('cond_median_spread'),
               frac_integer=float(np.mean(present == np.round(present))))
    if diag['frac_missing'] > MAX_MISSING:
        out['issues'].append(('warning', '%.0f%% of values missing' % (100 * diag['frac_missing'])))
    for a in actions:
        if 'transform' in a:
            out['issues'].append(('info', 'will be %s' % a))
    for w in warnings:
        out['issues'].append(('warning', w))
    return out


def validate(raw_dir, workers=None):
    """Return a Report for raw_dir."""
    report = Report()
    try:
        cfg = modconfig.load(raw_dir)
    except (modconfig.ConfigError, OSError, ValueError) as err:
        report.add('error', 'mod.toml', str(err))
        return report

    gene_path = os.path.join(raw_dir, 'gene_list.txt')
    if not os.path.exists(gene_path):
        report.add('error', 'gene_list.txt', 'missing')
        return report
    gene_ids, gene_names = read_gene_list(gene_path)
    if len(set(gene_ids)) < len(gene_ids):
        report.add('error', 'gene_list.txt', 'duplicate gene IDs')
    if len(set(gene_names)) < len(gene_names):
        report.add('error', 'gene_list.txt', 'duplicate gene names')
    name_to_index = {n: i for i, n in enumerate(gene_names)}

    ds_path = os.path.join(raw_dir, 'dataset_list.txt')
    if not os.path.exists(ds_path):
        report.add('error', 'dataset_list.txt', 'missing')
        return report
    entries = read_dataset_list(ds_path)
    report.next_dataset_id = max((i for i, _ in entries), default=0) + 1
    seen_ids, seen_files, unique = set(), {}, []
    for ds_id, fname in entries:
        if ds_id in seen_ids:
            report.add('error', 'dataset_list.txt', 'dataset ID %d used more than once' % ds_id)
        seen_ids.add(ds_id)
        if fname in seen_files:
            report.add('warning', fname, 'listed twice (IDs %d and %d); only ID %d is built'
                       % (seen_files[fname], ds_id, seen_files[fname]))
            continue
        seen_files[fname] = ds_id
        if not os.path.exists(os.path.join(raw_dir, 'expression', fname)):
            report.add('error', fname, 'listed in dataset_list.txt but not in expression/')
            continue
        unique.append((ds_id, fname))

    with ProcessPoolExecutor(workers) as pool:
        paths = [os.path.join(raw_dir, 'expression', f) for _, f in unique]
        results = pool.map(check_pcl, paths, [name_to_index] * len(paths), chunksize=8)
        for (ds_id, fname), res in zip(unique, results):
            for level, msg in res.pop('issues'):
                report.add(level, fname, msg)
            row = {'dataset_id': ds_id, 'filename': fname, **res}
            row['issues'] = '; '.join(i.message for i in report.issues if i.where == fname)
            report.datasets.append(row)

    _check_metadata(report, raw_dir, cfg, gene_names, [f for _, f in unique])
    return report


def _check_assays(report, cfg, meta, filenames):
    """Warn about missing or unknown assays and RNA-seq files that are not on a log scale."""
    where = cfg['datasets']['metadata'] or 'datasets'
    allowed = cfg['datasets']['assays']
    missing = [f for f in filenames if not (meta.get(f) or {}).get('assay')]
    if missing:
        report.add('warning', where, '%d datasets have no assay (they are left out when users filter by assay): '
                   '%s%s' % (len(missing), ', '.join(missing[:3]), ' …' if len(missing) > 3 else ''))
    unknown = sorted({m['assay'] for f in filenames if (m := meta.get(f)) and m.get('assay')} - set(allowed))
    for a in unknown:
        names = [f for f in filenames if (meta.get(f) or {}).get('assay') == a]
        report.add('warning', where, 'assay %r is not one of %s (%d datasets: %s%s)'
                   % (a, ', '.join(allowed), len(names), ', '.join(names[:3]), ' …' if len(names) > 3 else ''))
    stats = {d['filename']: d for d in report.datasets}
    no_type = []
    for f in filenames:
        m = meta.get(f) or {}
        if m.get('assay') != 'RNA-seq':
            continue
        if not m.get('value_type'):
            no_type.append(f)
        d = stats.get(f, {})
        if d.get('scale') != 'intensity':
            continue
        if d.get('frac_integer', 0) >= COUNT_LIKE:
            report.add('warning', f, 'RNA-seq values look like raw read counts; supply log2(TMM CPM + 1) '
                       'or log2(TPM + 1) instead (the build would only take log2(x + 1), with no '
                       'library-size normalisation)')
        else:
            report.add('warning', f, 'RNA-seq values are not log-scale (TPM, FPKM or CPM?); supply log2(TPM + 1) '
                       'or log2(TMM CPM + 1) (the build would take log2(x + 1) as supplied)')
    if no_type:
        report.add('warning', where, '%d RNA-seq datasets have no value_type (e.g. log2 TMM CPM, log2 TPM): %s%s'
                   % (len(no_type), ', '.join(no_type[:3]), ' …' if len(no_type) > 3 else ''))


def _check_metadata(report, raw_dir, cfg, gene_names, filenames):
    from spell_next.build import read_dataset_info, read_gene_info   # avoids a circular import

    genes = read_gene_info(raw_dir, cfg)
    if not genes:
        report.add('warning', cfg['genes']['file'] or 'genes', 'no gene names file: only systematic names searchable')
    else:
        unknown = sum(1 for n in gene_names if n not in genes)
        if unknown:
            report.add('warning', cfg['genes']['file'], '%d of %d genes in gene_list.txt have no entry'
                       % (unknown, len(gene_names)))
    meta = read_dataset_info(raw_dir, cfg)
    missing = [f for f in filenames if not (meta.get(f) or {}).get('title')]
    if missing:
        report.add('warning', cfg['datasets']['metadata'] or 'datasets',
                   '%d datasets have no metadata title (shown by file name): %s%s'
                   % (len(missing), ', '.join(missing[:5]), ' …' if len(missing) > 5 else ''))
    _check_assays(report, cfg, meta, filenames)
    extra = sorted(set(meta) - set(filenames))
    if extra:
        report.add('warning', cfg['datasets']['metadata'] or 'datasets',
                   '%d metadata rows name files not in dataset_list.txt: %s%s'
                   % (len(extra), ', '.join(extra[:5]), ' …' if len(extra) > 5 else ''))
    for g in [cfg['go']['obo']] + list(cfg['go']['gaf']):
        if not os.path.exists(modconfig.path(raw_dir, g)):
            report.add('warning', g, 'missing: GO enrichment will be disabled')
