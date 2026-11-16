"""Weekly reference update: spell update-refs fetches [updates] files and only swaps in good downloads."""
import gzip

import pytest

from spell_next import cli, modconfig


def run(*argv):
    return cli.main([str(a) for a in argv])


@pytest.fixture
def setup(tmp_path):
    """An instance whose mod.toml fetches two files from a local 'download site' (file:// URLs)."""
    root, site = tmp_path / 'instance', tmp_path / 'site'
    (root / 'raw' / 'meta').mkdir(parents=True)
    site.mkdir()
    (site / 'features.tab').write_text('S1\tORF\t\tYAL001C\tTFC3\n' * 50)
    with gzip.open(site / 'annot.gaf.gz', 'wt') as fh:
        fh.write('!gaf-version: 2.2\n' + 'SGD\tS1\tTFC3\t\tGO:0001\n' * 50)
    (root / 'raw' / 'mod.toml').write_text(
        '[updates]\n'
        '"meta/features.tab" = "%s"\n'
        '"meta/annot.gaf.gz" = "%s"\n' % ((site / 'features.tab').as_uri(), (site / 'annot.gaf.gz').as_uri()))
    return root, site


def test_first_fetch_then_unchanged_then_updated(setup, capsys):
    root, site = setup
    meta = root / 'raw' / 'meta'
    assert run('update-refs', root) == 0
    assert (meta / 'features.tab').read_text() == (site / 'features.tab').read_text()
    assert (meta / 'annot.gaf.gz').exists()

    assert run('update-refs', root) == 3                         # nothing changed: no rebuild needed
    assert 'All reference files are current' in capsys.readouterr().out

    old = (meta / 'features.tab').read_text()
    (site / 'features.tab').write_text(old.replace('TFC3', 'TFC3X'))
    assert run('update-refs', root) == 0
    assert 'TFC3X' in (meta / 'features.tab').read_text()
    assert (meta / 'features.tab.prev').read_text() == old       # previous version kept
    assert not list(meta.glob('*.download'))


@pytest.mark.parametrize('breakage, message', [
    (lambda site: (site / 'features.tab').write_text('S1\n'), 'under 50%'),
    (lambda site: (site / 'features.tab').write_text('<!DOCTYPE html><html>Not found</html>' * 100), 'HTML page'),
    (lambda site: (site / 'annot.gaf.gz').write_bytes((site / 'annot.gaf.gz').read_bytes()[:-20]),
     'not a complete gzip'),
    (lambda site: (site / 'features.tab').unlink(), 'FAILED'),
])
def test_bad_downloads_never_replace_current_files(setup, capsys, breakage, message):
    root, site = setup
    meta = root / 'raw' / 'meta'
    assert run('update-refs', root) == 0
    before = {p.name: p.read_bytes() for p in meta.iterdir()}
    breakage(site)
    assert run('update-refs', root) == 1
    assert message in capsys.readouterr().out
    assert {p.name: p.read_bytes() for p in meta.iterdir()} == before


def test_no_updates_section_is_a_no_op(tmp_path, capsys):
    (tmp_path / 'raw').mkdir()
    assert run('update-refs', tmp_path) == 3
    assert 'No [updates] section' in capsys.readouterr().out


@pytest.mark.parametrize('rel, url', [
    ('../outside.tab', 'https://example.org/x'),
    ('/etc/passwd', 'https://example.org/x'),
    ('meta/x.tab', 'ftp://example.org/x'),
])
def test_config_rejects_unsafe_updates(tmp_path, rel, url):
    (tmp_path / 'mod.toml').write_text('[updates]\n"%s" = "%s"\n' % (rel, url))
    with pytest.raises(modconfig.ConfigError):
        modconfig.load(str(tmp_path))
