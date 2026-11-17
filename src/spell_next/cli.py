"""Curator commands for one MOD instance: spell <command> ROOT ...

ROOT is the instance directory (or set SPELL_ROOT):

    ROOT/raw/       inputs curators edit: mod.toml, gene_list.txt, dataset_list.txt, expression/, meta/
    ROOT/builds/    one directory per build (YYYYMMDD-HHMMSS)
    ROOT/current    symlink to the build being served
    ROOT/logs/      monthly search logs (kept across builds)

Typical use: add PCL files to raw/expression/, list them in raw/dataset_list.txt and
raw/meta/datasets.tsv, then `spell validate ROOT` and `spell rebuild ROOT`. Gene names
and GO files listed under [updates] in mod.toml are refreshed by `spell update-refs ROOT`.
"""
import argparse
import io
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import time

from spell_next import modconfig

KEEP_BUILDS = 3


def paths(root):
    root = os.path.abspath(root)
    return {'root': root, 'raw': os.path.join(root, 'raw'), 'builds': os.path.join(root, 'builds'),
            'current': os.path.join(root, 'current'), 'logs': os.path.join(root, 'logs')}


def list_builds(p):
    if not os.path.isdir(p['builds']):
        return []
    return sorted(d for d in os.listdir(p['builds'])
                  if re.fullmatch(r'\d{8}-\d{6}', d) and os.path.isdir(os.path.join(p['builds'], d)))


def current_build(p):
    return os.path.basename(os.path.realpath(p['current'])) if os.path.islink(p['current']) else None


def switch(p, build):
    """Point ROOT/current at builds/BUILD atomically."""
    tmp = p['current'] + '.tmp'
    if os.path.lexists(tmp):
        os.remove(tmp)
    os.symlink(os.path.join('builds', build), tmp)
    os.replace(tmp, p['current'])


def restart(cmd):
    if cmd:
        print('restarting: %s' % cmd)
        subprocess.run(cmd, shell=True, check=True)
    elif os.environ.get('SPELL_ONEOFF') == '1':
        print('spellctl restarts the web app next.')
    else:
        print('Restart the web app to serve it (or pass --restart-cmd / set SPELL_RESTART_CMD).')


def build_summary(build_dir):
    db = sqlite3.connect('file:%s?mode=ro' % os.path.join(build_dir, 'spell.sqlite'), uri=True)
    try:
        info = dict(db.execute('SELECT key, value FROM build_info'))
        return {
            'built': info.get('built', ''),
            'datasets': db.execute('SELECT COUNT(*) FROM datasets').fetchone()[0],
            'genes': db.execute('SELECT COUNT(*) FROM genes').fetchone()[0],
            'conditions': db.execute('SELECT SUM(num_conds) FROM datasets').fetchone()[0],
            'files': {r[0] for r in db.execute('SELECT filename FROM datasets')},
        }
    finally:
        db.close()


def compare_builds(old_dir, new_dir, top=10):
    """How a new build differs from the live one: datasets added/removed, genes whose top hits changed."""
    old, new = build_summary(old_dir), build_summary(new_dir)
    db = sqlite3.connect('file:%s?mode=ro' % os.path.join(new_dir, 'spell.sqlite'), uri=True)
    db.execute("ATTACH DATABASE ? AS old", ('file:%s?mode=ro' % os.path.join(old_dir, 'spell.sqlite'),))
    changed = db.execute('''
        SELECT COUNT(*) FROM (
          SELECT n.gene_id, COUNT(o.hit_id) AS kept
          FROM (SELECT gene_id, hit_id FROM coexpression WHERE rank <= ?) n
          LEFT JOIN (SELECT gene_id, hit_id FROM old.coexpression WHERE rank <= ?) o
            ON o.gene_id = n.gene_id AND o.hit_id = n.hit_id
          GROUP BY n.gene_id) WHERE kept < ? - 2''', (top, top, top)).fetchone()[0]
    db.close()
    return {'added': sorted(new['files'] - old['files']), 'removed': sorted(old['files'] - new['files']),
            'genes_changed': changed, 'old': old, 'new': new}


def smoke_test(build_dir):
    """Fail early if a build is empty or cannot answer a search."""
    from spell_next.compendium import Compendium
    from spell_next.search import search

    s = build_summary(build_dir)
    if not s['datasets'] or not s['genes']:
        raise RuntimeError('build has %d datasets and %d genes' % (s['datasets'], s['genes']))
    db = sqlite3.connect(os.path.join(build_dir, 'spell.sqlite'))
    row = db.execute('SELECT gene_id, COUNT(*) FROM coexpression GROUP BY gene_id LIMIT 1').fetchone()
    db.close()
    comp = Compendium.load(os.path.join(build_dir, 'compendium.npz'))
    # Without precomputed results, test with a gene measured in the first dataset.
    gene = row[0] if row else comp.gene_ids[int(comp.datasets[0].gene_index[0])]
    result = search(comp, [gene], max_genes=5)
    if not result.dataset_ids:
        raise RuntimeError('a test search for gene %d found no datasets' % gene)


def print_report(report):
    for level in ('error', 'warning', 'info'):
        items = [i for i in report.issues if i.level == level]
        if items:
            print('\n%s (%d):' % (level.upper() + 'S', len(items)))
            for i in items[:200]:
                print('  %-34s %s' % (i.where[:34], i.message))
            if len(items) > 200:
                print('  … %d more in the report file' % (len(items) - 200))
    print('\n%d datasets checked: %d errors, %d warnings. Next free dataset ID: %d'
          % (len(report.datasets), len(report.errors), len(report.warnings), report.next_dataset_id))


def run_validate(p, workers):
    from spell_next.validate import validate

    report = validate(p['raw'], workers)
    print_report(report)
    out = os.path.join(p['root'], 'validation_report.tsv')
    report.write_tsv(out)
    print('Per-dataset report: %s' % out)
    return report


# ---- commands -------------------------------------------------------------------------

def cmd_init(args):
    p = paths(args.root)
    for d in (os.path.join(p['raw'], 'expression'), os.path.join(p['raw'], 'meta'), p['builds'], p['logs']):
        os.makedirs(d, exist_ok=True)
    target = os.path.join(p['raw'], 'mod.toml')
    if args.config and not os.path.exists(target):
        shutil.copy(args.config, target)
        print('copied %s to %s' % (args.config, target))
    print('Instance ready at %s. Next: fill raw/ (mod.toml, gene_list.txt, dataset_list.txt, expression/, '
          'meta/), then run: spell validate %s' % (p['root'], p['root']))


def cmd_validate(args):
    report = run_validate(paths(args.root), args.workers)
    return 1 if report.errors else 0


def cmd_rebuild(args):
    from spell_next.build import run_build

    p = paths(args.root)
    report = run_validate(p, args.workers)
    if report.errors:
        print('\nNot rebuilding: fix the errors above first.')
        return 1
    os.makedirs(p['builds'], exist_ok=True)
    name = time.strftime('%Y%m%d-%H%M%S')
    while os.path.exists(os.path.join(p['builds'], name)):      # never reuse a build name
        time.sleep(1)
        name = time.strftime('%Y%m%d-%H%M%S')
    partial = os.path.join(p['builds'], name + '.partial')
    try:
        run_build(p['raw'], partial, workers=args.workers, limit=args.limit)
        smoke_test(partial)
    except Exception:
        print('\nBuild failed; the live build is unchanged. Partial output: %s' % partial)
        raise
    shutil.copy(os.path.join(p['root'], 'validation_report.tsv'), partial)
    final = os.path.join(p['builds'], name)
    os.rename(partial, final)

    live = current_build(p)
    if live:
        diff = compare_builds(os.path.join(p['builds'], live), final)
        print('\nCompared with the live build %s: %d datasets added, %d removed; %d genes have more than 2 '
              'of their top-10 hits changed.' % (live, len(diff['added']), len(diff['removed']), diff['genes_changed']))
        for f in diff['added'][:20]:
            print('  + %s' % f)
        for f in diff['removed'][:20]:
            print('  - %s' % f)
    if args.no_switch:
        print('\nBuilt %s (not switched; use: spell rollback %s --to %s).' % (name, p['root'], name))
        return 0
    switch(p, name)
    print('\nNow serving build %s.' % name)
    prune(p, args.keep)
    restart(args.restart_cmd)
    return 0


def prune(p, keep):
    live = current_build(p)
    old = [b for b in list_builds(p) if b != live]
    for b in old[:max(0, len(old) - (keep - 1))]:
        shutil.rmtree(os.path.join(p['builds'], b))
        print('removed old build %s' % b)
    for d in os.listdir(p['builds']):
        if d.endswith('.partial'):
            shutil.rmtree(os.path.join(p['builds'], d), ignore_errors=True)


def cmd_rollback(args):
    p = paths(args.root)
    builds, live = list_builds(p), current_build(p)
    if args.to:
        if args.to not in builds:
            print('No build %s. Builds: %s' % (args.to, ', '.join(builds) or 'none'))
            return 1
        target = args.to
    else:
        older = [b for b in builds if live is None or b < live]
        if not older:
            print('No build older than %s to roll back to.' % live)
            return 1
        target = older[-1]
    switch(p, target)
    print('Now serving build %s (was %s).' % (target, live))
    restart(args.restart_cmd)
    return 0


def cmd_status(args):
    p = paths(args.root)
    live = current_build(p)
    builds = list_builds(p)
    try:
        cfg = modconfig.load(p['raw'])
        print('%s instance at %s' % (cfg['mod']['id'], p['root']))
    except Exception as err:  # noqa: BLE001 - status should still show builds
        print('instance at %s (config problem: %s)' % (p['root'], err))
    if not builds:
        print('No builds yet: run spell rebuild %s' % p['root'])
        return 0
    print('%-17s %-8s %9s %7s %11s  %s' % ('build', '', 'datasets', 'genes', 'conditions', 'built'))
    for b in reversed(builds):
        s = build_summary(os.path.join(p['builds'], b))
        print('%-17s %-8s %9d %7d %11d  %s' % (b, 'LIVE' if b == live else '', s['datasets'], s['genes'],
                                              s['conditions'] or 0, s['built']))
    return 0


def cmd_serve(args):
    import uvicorn

    p = paths(args.root)
    if not current_build(p):
        print('No live build: run spell rebuild %s' % p['root'])
        return 1
    os.environ['SPELL_NEXT_DATA'] = os.path.realpath(p['current'])
    os.environ.setdefault('SPELL_NEXT_LOG_DIR', p['logs'])
    uvicorn.run('spell_next.api:app', host=args.host, port=args.port)
    return 0


def cmd_migrate_metadata(args):
    """Turn the old SPELL metadata export into meta/datasets.tsv, which curators can edit."""
    from spell_next.build import read_dataset_info
    from spell_next.compendium import read_dataset_list

    p = paths(args.root)
    cfg = modconfig.load(p['raw'])
    if cfg['datasets']['format'] == 'dataset_table':
        print('Already using a dataset table: %s' % cfg['datasets']['metadata'])
        return 0
    meta = read_dataset_info(p['raw'], cfg)
    out_rel = 'meta/datasets.tsv'
    cols = ['filename', 'geo_id', 'pubmed_id', 'title', 'description', 'first_author', 'authors', 'journal',
            'year', 'tags', 'assay', 'value_type']
    rows = 0
    with open(os.path.join(p['raw'], out_rel), 'w') as out:
        out.write('\t'.join(cols) + '\n')
        listed = dict.fromkeys(f for _, f in read_dataset_list(os.path.join(p['raw'], 'dataset_list.txt')))
        for fname in list(listed) + sorted(set(meta) - set(listed)):
            m = dict(meta.get(fname, {}))
            m['filename'] = fname
            geo = re.match(r'(GSE\d+)', fname)
            m['geo_id'] = m.get('geo_id') or (geo.group(1) if geo else '')
            out.write('\t'.join(re.sub(r'[\t\r\n]+', ' ', str(m.get(c) or '')) for c in cols) + '\n')
            rows += 1
    print('wrote %d rows to raw/%s' % (rows, out_rel))
    toml_path = os.path.join(p['raw'], 'mod.toml')
    text = open(toml_path).read() if os.path.exists(toml_path) else ''
    section = re.search(r'(?ms)^\[datasets\]\n(.*?)(?=^\[|\Z)', text)
    new_section = '[datasets]\nmetadata = "%s"\nformat = "dataset_table"\n\n' % out_rel
    text = text[:section.start()] + new_section + text[section.end():] if section else text + '\n' + new_section
    with open(toml_path, 'w') as fh:
        fh.write(text)
    print('updated %s: [datasets] now reads raw/%s' % (toml_path, out_rel))
    return 0


def _dataset_geo(p, cfg):
    """{filename: GEO series} from the dataset table, or from file names where it has none."""
    from spell_next.build import read_dataset_info
    from spell_next.geo import geo_id_of

    meta = read_dataset_info(p['raw'], cfg)
    return lambda f: (meta.get(f) or {}).get('geo_id') or geo_id_of(f)


def _fingerprint_inputs(p, listed, retired, workers):
    """Fingerprints of listed and retired datasets: (name_to_index, [(label, path, fp)])."""
    from spell_next import duplicates
    from spell_next.compendium import read_gene_list

    _, names = read_gene_list(os.path.join(p['raw'], 'gene_list.txt'))
    name_to_index = {n: i for i, n in enumerate(names)}
    items = []
    for folder, entries, tag in ((os.path.join(p['raw'], 'expression'), listed, ''),
                                 (os.path.join(p['raw'], 'retired', 'expression'), retired, 'retired ')):
        for ds_id, f in entries:
            if os.path.exists(os.path.join(folder, f)):
                items.append(('%sdataset %d (%s)' % (tag, ds_id, f), os.path.join(folder, f)))
    fps = duplicates.fingerprints(p['raw'], [path for _, path in items], name_to_index, workers)
    return name_to_index, [(label, path, fps[path]) for label, path in items]


def _report_duplicates(p, cfg, files, gse, listed, retired, workers):
    """Print possible duplicates of the new files; return the new files whose data is already in SPELL."""
    from spell_next import duplicates

    print('\nChecking for duplicates (first run reads every dataset; later runs only new or changed files)...')
    name_to_index, existing = _fingerprint_inputs(p, listed, retired, workers)
    paths = [os.path.join(p['raw'], 'expression', f) for f in files]
    fps = duplicates.fingerprints(p['raw'], paths, name_to_index, workers)
    geo_existing = _dataset_geo(p, cfg)

    def geo_of(label):
        if label in gse:
            return gse[label]
        return geo_existing(label.rsplit('(', 1)[-1].rstrip(')'))

    found = duplicates.compare([(f, path, fps[path]) for f, path in zip(files, paths)], existing, geo_of)
    same_data = []
    for f, hits in found.items():
        for kind, other in hits:
            serious = kind == 'same data' and not other.startswith('retired ')
            print('  %s %s: %s as %s' % ('DUPLICATE' if serious else 'note:', f, kind, other))
        if any(k == 'same data' and not o.startswith('retired ') for k, o in hits):
            same_data.append(f)
    if not found:
        print('  none found')
    return same_data


def cmd_find_duplicates(args):
    """Report datasets in SPELL that share data, condition names or a GEO series."""
    from spell_next import duplicates
    from spell_next.compendium import read_dataset_list

    p = paths(args.root)
    cfg = modconfig.load(p['raw'])
    listed = read_dataset_list(os.path.join(p['raw'], 'dataset_list.txt'))
    seen, unique = set(), []
    for ds_id, f in listed:                   # a file listed twice is reported by validate already
        if f not in seen:
            seen.add(f)
            unique.append((ds_id, f))
    print('Fingerprinting %d datasets (cached after the first run)...' % len(unique))
    _, items = _fingerprint_inputs(p, unique, [], args.workers)
    geo_of = _dataset_geo(p, cfg)
    found = duplicates.groups(items, lambda label: geo_of(label.rsplit('(', 1)[-1].rstrip(')')))
    reported = set()
    for kind, labels in found:
        key = frozenset(labels)
        if key in reported:                   # e.g. same data implies same conditions: report once
            continue
        reported.add(key)
        print('\n%s (%d datasets):' % (kind, len(labels)))
        for label in labels:
            print('  ' + label)
    if not reported:
        print('No duplicates found.')
    print('\n"same data" means identical values; "same conditions" the same column names (often the same '
          'experiment processed differently); "same GEO series" can be legitimate (one series split into '
          'several datasets).')
    return 0


def _assign_ids(files, supplied, listed, retired):
    """Dataset IDs for new files: the table's dataset_id where given (e.g. to match the old SPELL), else the
    next free IDs. IDs of listed or retired datasets are never reused. Returns ({file: id}, problems)."""
    used = {i: f for i, f in retired}
    used.update({i: f for i, f in listed})
    ids, problems = {}, []
    for f in files:
        raw = (supplied.get(f) or {}).get('dataset_id') or (supplied.get(f) or {}).get('id') or ''
        if not raw:
            continue
        if not raw.isdigit() or int(raw) == 0:
            problems.append('%s: dataset_id %r is not a positive whole number' % (f, raw))
        elif int(raw) in used:
            problems.append('%s: dataset_id %s is already used by %s%s' % (
                f, raw, used[int(raw)], ' (retired)' if int(raw) in dict(retired) and
                int(raw) not in dict(listed) else ''))
        elif int(raw) in ids.values():
            problems.append('%s: dataset_id %s is given to two files' % (f, raw))
        else:
            ids[f] = int(raw)
    taken = set(used) | set(ids.values())
    nxt = max(taken, default=0) + 1
    for f in files:
        if f not in ids:
            ids[f] = nxt
            nxt += 1
    return ids, problems


FETCHED = ('title', 'pubmed_id', 'description', 'first_author', 'authors', 'journal', 'year', 'assay')


def read_metadata_table(path):
    """A curator's metadata table (tab- or comma-separated, header row with a filename column):
    {filename: {column: value}}, columns lower-cased."""
    import csv

    with open(path, newline='', encoding='utf-8-sig', errors='replace') as fh:
        text = fh.read()
    delimiter = '\t' if '\t' in text.split('\n', 1)[0] else ','
    out = {}
    for row in csv.DictReader(io.StringIO(text), delimiter=delimiter):
        row = {(k or '').strip().lower(): re.sub(r'[\t\r\n]+', ' ', (v or '')).strip() for k, v in row.items()}
        name = os.path.basename(row.get('filename') or '')
        if name:
            out[name] = row
    if not out:
        raise ValueError('no rows with a filename column')
    return out


def _backup(path):
    if os.path.exists(path):
        shutil.copy2(path, '%s.bak-%s' % (path, time.strftime('%Y%m%d-%H%M%S')))


def cmd_add_datasets(args):
    """List new expression files in dataset_list.txt and datasets.tsv. Metadata comes from --metadata (a
    curator's table) where given, and from GEO and PubMed for whatever that leaves blank."""
    import csv

    from spell_next import geo
    from spell_next.compendium import read_dataset_list

    p = paths(args.root)
    cfg = modconfig.load(p['raw'])
    if cfg['datasets']['format'] != 'dataset_table' or not cfg['datasets']['metadata']:
        print('add-datasets needs a dataset table ([datasets] format = "dataset_table"); '
              'run spell migrate-metadata first.')
        return 1
    list_path = os.path.join(p['raw'], 'dataset_list.txt')
    meta_path = os.path.join(p['raw'], cfg['datasets']['metadata'])
    expr = os.path.join(p['raw'], 'expression')
    listed = read_dataset_list(list_path) if os.path.exists(list_path) else []
    retired_path = os.path.join(p['raw'], 'retired', 'dataset_list.txt')
    retired = read_dataset_list(retired_path) if os.path.exists(retired_path) else []
    listed_files = {f for _, f in listed}

    supplied = {}
    if args.metadata:
        try:
            supplied = read_metadata_table(args.metadata)
        except (OSError, ValueError) as err:
            print('Cannot read %s: %s' % (args.metadata, err))
            return 1
    if args.files:
        files = [os.path.basename(f) for f in args.files]
    elif supplied:                              # the files the metadata table describes
        files = list(supplied)
    else:                                       # every expression file not yet listed
        files = sorted(f for f in os.listdir(expr) if not f.endswith('.sinfo') and not f.startswith('.')
                       and f not in listed_files)
    problems = [f + ': not in raw/expression/' for f in files if not os.path.exists(os.path.join(expr, f))]
    problems += [f + ': already in dataset_list.txt' for f in files if f in listed_files]
    if problems:
        print('Nothing added:\n  ' + '\n  '.join(problems))
        return 1
    if not files:
        print('No new files in raw/expression/: every file there is already listed.')
        return 0

    # Display-only datasets have IDs too (in the dataset table); never hand those out again.
    from spell_next.build import read_dataset_info
    display = [(m['dataset_id'], f) for f, m in read_dataset_info(p['raw'], cfg).items()
               if m.get('display_only') and m.get('dataset_id') is not None]
    ids, problems = _assign_ids(files, supplied, listed, retired + display)
    if problems:
        print('Nothing added:\n  ' + '\n  '.join(problems))
        return 1

    gse = {f: (supplied.get(f) or {}).get('geo_id') or geo.geo_id_of(f) for f in files}
    # Look up only datasets the curator's table leaves incomplete.
    lookup = {g for f, g in gse.items() if g and any(not (supplied.get(f) or {}).get(c) for c in FETCHED)}
    meta = {}
    if not args.no_fetch and lookup:
        print('Fetching metadata for %d GEO series from NCBI...' % len(lookup))
        try:
            meta = geo.metadata_for(sorted(lookup))
        except geo.NCBIError as err:
            print('NCBI lookup failed (%s); nothing added. Retry, or use --no-fetch to add rows without '
                  'metadata.' % err)
            return 1

    with open(meta_path, newline='', encoding='utf-8', errors='replace') as fh:
        cols = next(csv.reader(fh, delimiter='\t'), None) or []
    header = list(cols)
    for c in ('filename', 'geo_id', 'pubmed_id', 'title', 'description', 'first_author', 'authors', 'journal',
              'year', 'tags', 'assay', 'value_type'):
        if c not in cols:
            cols.append(c)
    rows, lines, notes = [], [], []
    for f in files:
        m = dict(meta.get(gse[f]) or {})
        found = bool(m) or gse[f] not in lookup
        m.pop('gdstype', None)
        m.update(filename=f, geo_id=gse[f] or '', tags=args.tags or '',
                 assay=args.assay or m.get('assay', ''), value_type=args.value_type or '')
        m.update({c: v for c, v in (supplied.get(f) or {}).items() if v})    # the curator's values win
        rows.append({c: re.sub(r'[\t\r\n]+', ' ', str(m.get(c) or '')) for c in cols})
        lines.append('%d\t%s\n' % (ids[f], f))
        if not gse[f]:
            notes.append('%s: no GSE ID in the file name, so no metadata was fetched' % f)
        elif not found and not args.no_fetch:
            notes.append('%s: %s not found in GEO' % (f, gse[f]))
        elif not m.get('pubmed_id'):
            notes.append('%s: GEO lists no PubMed record; title is the GEO series title' % f)
        if not rows[-1]['assay']:
            notes.append('%s: assay unknown; set it in %s' % (f, cfg['datasets']['metadata']))

    print('\n%-6s %-34s %-11s %-10s %-11s %s' % ('ID', 'file', 'GEO', 'PubMed', 'assay', 'title'))
    for line, r in zip(lines, rows):
        print('%-6s %-34s %-11s %-10s %-11s %s' % (line.split('\t')[0], r['filename'][:34], r['geo_id'],
                                                   r['pubmed_id'], r['assay'], r['title'][:60]))
    for n in notes:
        print('  note: ' + n)
    same_data = _report_duplicates(p, cfg, files, gse, listed, retired, args.workers)
    if args.dry_run:
        print('\nDry run: nothing written.')
        return 0
    if same_data and not args.allow_duplicates:
        print('\nNothing added: %d file(s) hold the same data as a dataset already in SPELL (or earlier in this '
              'batch). Remove them, or rerun with --allow-duplicates to add them anyway.' % len(same_data))
        return 1
    _backup(list_path)
    _backup(meta_path)
    with open(list_path, 'a') as fh:
        fh.writelines(lines)
    if header != cols[:len(header)] or len(header) < len(cols):
        # New columns: rewrite the table with the full header, keeping every existing row.
        with open(meta_path, newline='', encoding='utf-8', errors='replace') as fh:
            old = list(csv.DictReader(fh, delimiter='\t'))
        with open(meta_path, 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, cols, delimiter='\t', lineterminator='\n', restval='')
            w.writeheader()
            w.writerows(old + rows)
    else:
        with open(meta_path, 'rb') as fh:
            fh.seek(-1, os.SEEK_END)
            missing_newline = fh.read(1) != b'\n'
        with open(meta_path, 'a', newline='', encoding='utf-8') as fh:
            if missing_newline:
                fh.write('\n')
            csv.DictWriter(fh, cols, delimiter='\t', lineterminator='\n').writerows(rows)
    print('\nAdded %d datasets (IDs %d-%d) to dataset_list.txt and %s (backups kept as *.bak-*).'
          % (len(files), min(ids.values()), max(ids.values()), cfg['datasets']['metadata']))
    print('Next: check the new rows (add tags; set value_type for RNA-seq), then validate and rebuild.')
    return 0


def cmd_update_refs(args):
    """Fetch the reference files listed under [updates]; exit 0 if any changed, 3 if none did, 1 on failure."""
    from spell_next import refs

    p = paths(args.root)
    cfg = modconfig.load(p['raw'])
    if not cfg['updates']:
        print('No [updates] section in raw/mod.toml: nothing to fetch.')
        return 3
    print('Fetching %d reference file(s):' % len(cfg['updates']))
    changed, failed = refs.update(p['raw'], cfg['updates'])
    if failed:
        print('%d download(s) failed; those files were left as they were.' % len(failed))
        return 1
    if not changed:
        print('All reference files are current.')
        return 3
    print('%d file(s) updated; rebuild to use them.' % len(changed))
    return 0


def main(argv=None):
    # Workers are processes; keep BLAS single-threaded in each (set before numpy is imported).
    for var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(var, '1')
    ap = argparse.ArgumentParser(prog='spell', description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='command', required=True)
    root_default = os.environ.get('SPELL_ROOT')

    def add(name, func, help_text):
        sp = sub.add_parser(name, help=help_text, description=help_text)
        sp.add_argument('root', nargs='?' if root_default else None, default=root_default,
                        help='instance directory (default: $SPELL_ROOT)')
        sp.set_defaults(func=func)
        return sp

    sp = add('init', cmd_init, 'create the directory layout for a new instance')
    sp.add_argument('--config', help='config template to copy to raw/mod.toml (e.g. configs/wb.example.toml)')
    sp = add('validate', cmd_validate, 'check raw inputs without building')
    sp.add_argument('--workers', type=int)
    sp = add('rebuild', cmd_rebuild, 'validate, build, check and switch to the new build')
    sp.add_argument('--workers', type=int)
    sp.add_argument('--limit', type=int, help=argparse.SUPPRESS)
    sp.add_argument('--keep', type=int, default=KEEP_BUILDS, help='builds to keep, live one included')
    sp.add_argument('--no-switch', action='store_true', help='build but keep serving the current build')
    sp.add_argument('--restart-cmd', default=os.environ.get('SPELL_RESTART_CMD'))
    sp = add('rollback', cmd_rollback, 'serve the previous build again (or --to BUILD)')
    sp.add_argument('--to')
    sp.add_argument('--restart-cmd', default=os.environ.get('SPELL_RESTART_CMD'))
    add('status', cmd_status, 'show the live build and the build history')
    sp = add('serve', cmd_serve, 'run the web app on the live build')
    sp.add_argument('--host', default='0.0.0.0')
    sp.add_argument('--port', type=int, default=8000)
    add('migrate-metadata', cmd_migrate_metadata, 'convert the old SPELL metadata export into meta/datasets.tsv')
    add('update-refs', cmd_update_refs, 'download the reference files listed under [updates] in raw/mod.toml')
    sp = add('add-datasets', cmd_add_datasets,
             'add new expression files to the dataset lists, with metadata from GEO and PubMed')
    sp.add_argument('files', nargs='*', help='files in raw/expression/ (default: those in --metadata, else '
                                             'every file not yet listed)')
    sp.add_argument('--metadata', help='table of the new datasets (filename plus any datasets.tsv columns, and '
                                       'optionally dataset_id); its values are used as given, GEO/PubMed fill '
                                       'blanks')
    sp.add_argument('--assay', help='assay for all of them (default: from GEO, e.g. microarray or RNA-seq)')
    sp.add_argument('--value-type', help='value_type for all of them, e.g. "log2 TMM CPM"')
    sp.add_argument('--tags', help='tags for all of them, |-separated')
    sp.add_argument('--dry-run', action='store_true', help='show what would be added, write nothing')
    sp.add_argument('--no-fetch', action='store_true', help='do not look up GEO/PubMed (rows get file names only)')
    sp.add_argument('--allow-duplicates', action='store_true',
                    help='add files even if their data is already in SPELL')
    sp.add_argument('--workers', type=int)
    sp = add('find-duplicates', cmd_find_duplicates,
             'list datasets in SPELL that share data, condition names or a GEO series')
    sp.add_argument('--workers', type=int)

    args = ap.parse_args(argv)
    return args.func(args) or 0


if __name__ == '__main__':
    sys.exit(main())
