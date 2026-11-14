"""Curator commands for one MOD instance: spell <command> ROOT ...

ROOT is the instance directory (or set SPELL_ROOT):

    ROOT/raw/       inputs curators edit: mod.toml, gene_list.txt, dataset_list.txt, expression/, meta/
    ROOT/builds/    one directory per build (YYYYMMDD-HHMMSS)
    ROOT/current    symlink to the build being served
    ROOT/logs/      monthly search logs (kept across builds)

Typical use: add PCL files to raw/expression/, list them in raw/dataset_list.txt and
raw/meta/datasets.tsv, then `spell validate ROOT` and `spell rebuild ROOT`.
"""
import argparse
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
    gene, n = db.execute('SELECT gene_id, COUNT(*) FROM coexpression GROUP BY gene_id LIMIT 1').fetchone()
    db.close()
    comp = Compendium.load(os.path.join(build_dir, 'compendium.npz'))
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
            'year', 'tags']
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

    args = ap.parse_args(argv)
    return args.func(args) or 0


if __name__ == '__main__':
    sys.exit(main())
