"""spell add-datasets: new expression files get IDs and metadata rows, filled from GEO and PubMed (faked here)."""
import csv

import pytest

from conftest import make_toy_raw
from spell_next import cli, geo

PUBMED_XML = (
    '<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>111</PMID><Article><Abstract>'
    '<AbstractText Label="BACKGROUND">Cells grow.</AbstractText>'
    '<AbstractText Label="RESULTS">They <i>divide</i>.</AbstractText>'
    '</Abstract></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>')
CALLS = []


def fake_ncbi(tool, **params):
    """Stand-in for geo._post: GSE1001 is RNA-seq with a paper, GSE1002 an array with none."""
    CALLS.append((tool, params.get('db'), params.get('term') or params.get('id')))
    if tool == 'esearch.fcgi':
        return {'esearchresult': {'idlist': ['2001', '2002']}}
    if tool == 'esummary.fcgi' and params['db'] == 'gds':
        return {'result': {'uids': ['2001', '2002'],
                           '2001': {'accession': 'GSE1001', 'title': 'Series one', 'summary': 'GEO summary one',
                                    'pubmedids': ['111'],
                                    'gdstype': 'Expression profiling by high throughput sequencing'},
                           '2002': {'accession': 'GSE1002', 'title': 'Series two', 'summary': 'GEO summary two',
                                    'pubmedids': [], 'gdstype': 'Expression profiling by array'}}}
    if tool == 'esummary.fcgi' and params['db'] == 'pubmed':
        return {'result': {'uids': ['111'], '111': {
            'title': 'Paper one.', 'pubdate': '2025 Mar 3', 'fulljournalname': 'Yeast',
            'authors': [{'name': 'Doe J', 'authtype': 'Author'}, {'name': 'Roe R', 'authtype': 'Author'},
                        {'name': 'Yeast Consortium', 'authtype': 'CollectiveName'}]}}}
    if tool == 'efetch.fcgi':
        return PUBMED_XML
    raise AssertionError(tool)


@pytest.fixture
def root(tmp_path, monkeypatch):
    root = tmp_path / 'instance'
    assert cli.main(['init', str(root)]) == 0
    make_toy_raw(root / 'raw')
    pcl = (root / 'raw' / 'expression' / 'A_set0.pcl').read_text()
    for name in ('GSE1001.final.pcl', 'GSE1002_set0.pcl', 'local.pcl'):
        (root / 'raw' / 'expression' / name).write_text(pcl)
    CALLS.clear()
    monkeypatch.setattr(geo, '_post', fake_ncbi)
    monkeypatch.setattr(geo.metadata_for, '__defaults__', (fake_ncbi,))
    return root


def rows(root):
    with open(root / 'raw' / 'meta' / 'datasets.tsv', newline='') as fh:
        return {r['filename']: r for r in csv.DictReader(fh, delimiter='\t')}


def test_adds_unlisted_files_with_metadata(root, capsys):
    assert cli.main(['add-datasets', str(root)]) == 0
    listing = (root / 'raw' / 'dataset_list.txt').read_text().splitlines()
    assert listing[-3:] == ['5\tGSE1001.final.pcl', '6\tGSE1002_set0.pcl', '7\tlocal.pcl']
    r = rows(root)
    one = r['GSE1001.final.pcl']
    assert (one['pubmed_id'], one['title'], one['first_author'], one['authors'], one['journal'], one['year']) == \
        ('111', 'Paper one.', 'Doe J', 'Doe J, Roe R', 'Yeast', '2025')
    assert one['description'] == 'BACKGROUND: Cells grow. RESULTS: They divide.'   # the abstract, not GEO's
    assert one['assay'] == 'RNA-seq'
    two = r['GSE1002_set0.pcl']
    assert (two['title'], two['description'], two['assay'], two['pubmed_id']) == \
        ('Series two', 'GEO summary two', 'microarray', '')
    assert r['local.pcl']['title'] == '' and r['local.pcl']['geo_id'] == ''
    out = capsys.readouterr().out
    assert 'no GSE ID in the file name' in out and 'GEO lists no PubMed record' in out
    assert list((root / 'raw').glob('dataset_list.txt.bak-*')) and list((root / 'raw' / 'meta').glob('*.bak-*'))
    assert r['A_set0.pcl']['title'] == 'Alpha worm study 0'              # existing rows kept intact
    assert cli.main(['validate', str(root), '--workers', '2']) == 0         # still a valid instance


def test_options_dry_run_and_refusals(root, capsys):
    before = (root / 'raw' / 'dataset_list.txt').read_text()
    assert cli.main(['add-datasets', str(root), 'GSE1001.final.pcl', '--dry-run']) == 0
    assert (root / 'raw' / 'dataset_list.txt').read_text() == before
    assert cli.main(['add-datasets', str(root), 'GSE1001.final.pcl', '--assay', 'microarray',
                     '--value-type', 'log2 TPM', '--tags', 'stress|heat shock']) == 0
    one = rows(root)['GSE1001.final.pcl']
    assert (one['assay'], one['value_type'], one['tags']) == ('microarray', 'log2 TPM', 'stress|heat shock')
    assert cli.main(['add-datasets', str(root), 'GSE1001.final.pcl']) == 1      # already listed
    assert cli.main(['add-datasets', str(root), 'missing.pcl']) == 1            # not in expression/
    assert 'already in dataset_list.txt' in capsys.readouterr().out


def test_curators_metadata_table_is_used_as_given(root, tmp_path):
    """Edith's workflow: her table has full details for GSE1001 and only tags/assay for GSE1002."""
    table = tmp_path / 'new_datasets.csv'
    table.write_text(
        'filename,geo_id,pubmed_id,title,first_author,authors,journal,year,description,tags,assay,value_type\n'
        'GSE1001.final.pcl,GSE1001,999,My title,Smith A,"Smith A, Jones B",Genetics,2024,My summary,'
        'stress,RNA-seq,log2 TMM CPM\n'
        'GSE1002_set0.pcl,,,,,,,,,growth,microarray,log2 ratio\n')
    assert cli.main(['add-datasets', str(root), '--metadata', str(table)]) == 0
    r = rows(root)
    assert set(r) >= {'GSE1001.final.pcl', 'GSE1002_set0.pcl'} and 'local.pcl' not in r   # only the table's files
    one = r['GSE1001.final.pcl']
    assert (one['pubmed_id'], one['title'], one['authors'], one['journal'], one['description'], one['value_type']) \
        == ('999', 'My title', 'Smith A, Jones B', 'Genetics', 'My summary', 'log2 TMM CPM')
    two = r['GSE1002_set0.pcl']
    assert (two['title'], two['tags'], two['assay'], two['value_type']) == \
        ('Series two', 'growth', 'microarray', 'log2 ratio')           # GEO fills the blanks
    searched = ' '.join(str(c[2]) for c in CALLS if c[0] == 'esearch.fcgi')
    assert 'GSE1002' in searched and 'GSE1001' not in searched        # complete rows are not looked up


def test_ids_skip_retired_datasets(root):
    (root / 'raw' / 'retired').mkdir()
    (root / 'raw' / 'retired' / 'dataset_list.txt').write_text('40\told.pcl\n')
    assert cli.main(['add-datasets', str(root), 'local.pcl', '--no-fetch']) == 0
    assert (root / 'raw' / 'dataset_list.txt').read_text().splitlines()[-1] == '41\tlocal.pcl'


def test_assay_from_geo_type():
    assert geo.assay_of('Expression profiling by high throughput sequencing') == 'RNA-seq'
    assert geo.assay_of('Expression profiling by array; Other') == 'microarray'
    assert geo.assay_of('Other') is None
    assert geo.geo_id_of('GSE28213GPL13340.sfp.pcl') == 'GSE28213' and geo.geo_id_of('2010.Gasch00.pcl') is None
