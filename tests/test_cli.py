"""Curator workflow: init, validate, rebuild, add a dataset, rollback, prune, migrate metadata."""
import os

import numpy as np
import pytest

from conftest import N_CONDS, make_toy_raw
from spell_next import cli, modconfig


def run(*argv):
    return cli.main([str(a) for a in argv])


def live(root):
    return cli.current_build(cli.paths(root))


def add_dataset(raw, name, ds_id, title, n_conds=N_CONDS, rows=None):
    rng = np.random.default_rng(ds_id)
    lines = ['ID\tNAME\tGWEIGHT\t' + '\t'.join('c%d' % c for c in range(n_conds))]
    for i in range(1, 13):
        lines.append('A%d\tA%d\t1\t' % (i, i) + '\t'.join('%.3f' % v for v in 8 + rng.normal(size=n_conds)))
    if rows:
        lines.extend(rows)
    (raw / 'expression' / name).write_text('\n'.join(lines) + '\n')
    with open(raw / 'dataset_list.txt', 'a') as fh:
        fh.write('%d\t%s\n' % (ds_id, name))
    with open(raw / 'meta' / 'datasets.tsv', 'a') as fh:
        fh.write('%s\tGSE%d\t%s\tCurator C\t2026\tstress\n' % (name, 8000 + ds_id, title))


@pytest.fixture
def root(tmp_path):
    root = tmp_path / 'instance'
    assert run('init', root) == 0
    make_toy_raw(root / 'raw')
    return root


def test_init_creates_layout_and_copies_config(tmp_path):
    root = tmp_path / 'wb'
    template = os.path.join(os.path.dirname(__file__), '..', 'configs', 'wb.example.toml')
    assert run('init', root, '--config', template) == 0
    for d in ('raw/expression', 'raw/meta', 'builds', 'logs'):
        assert (root / d).is_dir()
    assert 'WormBase' in (root / 'raw' / 'mod.toml').read_text()


def test_validate_passes_on_good_inputs(root, capsys):
    assert run('validate', root, '--workers', '2') == 0
    out = capsys.readouterr().out
    assert '4 datasets checked: 0 errors' in out and 'Next free dataset ID: 5' in out
    assert (root / 'validation_report.tsv').exists()


@pytest.mark.parametrize('breakage, message', [
    ('duplicate_id', 'dataset ID 1 used more than once'),
    ('missing_file', 'not in expression/'),
    ('too_many_values', 'cannot be read'),
    ('unknown_genes', 'no rows match a gene'),
])
def test_validate_reports_errors(root, capsys, breakage, message):
    raw = root / 'raw'
    if breakage == 'duplicate_id':
        with open(raw / 'dataset_list.txt', 'a') as fh:
            fh.write('1\tA_set1.pcl\n')
    elif breakage == 'missing_file':
        with open(raw / 'dataset_list.txt', 'a') as fh:
            fh.write('9\tnot_there.pcl\n')
    elif breakage == 'too_many_values':
        add_dataset(raw, 'bad.pcl', 9, 'Bad', rows=['A1\tA1\t1\t' + '\t'.join(['1'] * (N_CONDS + 3))])
    else:
        (raw / 'expression' / 'alien.pcl').write_text('ID\tNAME\tc1\tc2\nX1\tX1\t1\t2\nX2\tX2\t3\t4\n')
        with open(raw / 'dataset_list.txt', 'a') as fh:
            fh.write('9\talien.pcl\n')
    assert run('validate', root, '--workers', '2') == 1
    assert message in capsys.readouterr().out


def test_rebuild_add_dataset_rollback_and_prune(root, capsys):
    restart_marker = root / 'restarted'
    assert run('rebuild', root, '--workers', '2', '--restart-cmd', 'touch %s' % restart_marker) == 0
    first = live(root)
    assert first and restart_marker.exists()
    assert os.path.exists(root / 'current' / 'spell.sqlite')
    assert (root / 'builds' / first / 'validation_report.tsv').exists()

    add_dataset(root / 'raw', 'A_new.pcl', 5, 'A brand new study')
    assert run('rebuild', root, '--workers', '2') == 0
    second = live(root)
    assert second > first
    out = capsys.readouterr().out
    assert '1 datasets added, 0 removed' in out and '+ A_new.pcl' in out
    assert cli.build_summary(str(root / 'current'))['datasets'] == 5

    assert run('rollback', root) == 0
    assert live(root) == first
    assert run('rollback', root, '--to', second) == 0
    assert live(root) == second

    assert run('rebuild', root, '--workers', '2', '--keep', '2') == 0
    third = live(root)
    assert cli.list_builds(cli.paths(root)) == [second, third]      # oldest pruned

    assert run('status', root) == 0
    assert 'LIVE' in capsys.readouterr().out


def test_failed_validation_leaves_live_build_alone(root, capsys):
    assert run('rebuild', root, '--workers', '2') == 0
    before = live(root)
    with open(root / 'raw' / 'dataset_list.txt', 'a') as fh:
        fh.write('9\tnot_there.pcl\n')
    assert run('rebuild', root, '--workers', '2') == 1
    assert live(root) == before
    assert 'Not rebuilding' in capsys.readouterr().out


def test_no_switch_builds_without_serving(root):
    assert run('rebuild', root, '--workers', '2') == 0
    before = live(root)
    assert run('rebuild', root, '--workers', '2', '--no-switch') == 0
    assert live(root) == before and len(cli.list_builds(cli.paths(root))) == 2


def test_new_datasets_need_no_sinfo_file(root):
    assert not list((root / 'raw' / 'expression').glob('A_new*'))
    add_dataset(root / 'raw', 'A_new.pcl', 5, 'No sinfo')
    assert not (root / 'raw' / 'expression' / 'A_new.pcl.sinfo').exists()
    assert run('rebuild', root, '--workers', '2') == 0


def test_migrate_metadata_from_old_spell_export(root):
    raw = root / 'raw'
    (raw / 'meta' / 'SQLdump2023.txt').write_text(
        '123\tA_set0.pcl\tOld title\tOld summary\t10\tSmith J\tSmith J, Doe A\tOld title\tJ Test\t2010\tstress\t1\n')
    toml = (raw / 'mod.toml').read_text().replace(
        '[datasets]\nmetadata = "meta/datasets.tsv"\nformat = "dataset_table"',
        '[datasets]\nmetadata = "meta/SQLdump2023.txt"\nformat = "spell_dump"')
    (raw / 'mod.toml').write_text(toml)
    assert run('migrate-metadata', root) == 0
    cfg = modconfig.load(str(raw))
    assert cfg['datasets'] == {'metadata': 'meta/datasets.tsv', 'format': 'dataset_table'}
    rows = (raw / 'meta' / 'datasets.tsv').read_text().splitlines()
    assert rows[0].startswith('filename\tgeo_id\tpubmed_id\ttitle')
    assert any(r.startswith('A_set0.pcl\t\t123\tOld title\tOld summary\tSmith J') for r in rows)
    assert len(rows) == 5                       # header + the 4 listed datasets
    assert run('validate', root, '--workers', '2') == 0
