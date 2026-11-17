"""Display-only datasets (metadata and links, no expression data) and dataset links."""
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import make_toy_raw
from spell_next import validate

DISPLAY = ('1001\tSC1.sc.paper\t\t\tA single-cell atlas\tAtlas of cells.\tCell A\t\t\t2025\tstress\t'
           'Single-Cell RNA-Seq\t\tAlpha worm\tyes\tGEO Record GSE77 <https://example.org/geo/GSE77> | '
           'Lab Atlas <https://atlas.example.org/>\n')


def add_display_rows(raw):
    tsv = raw / 'meta' / 'datasets.tsv'
    lines = tsv.read_text().splitlines()
    header = lines[0] + '\tdataset_id\tspecies\tdisplay_only\tlinks'
    rows = [line + '\t\t\t\t' for line in lines[1:]]
    rows[0] = rows[0][:-4] + '\t\t\t\tPaper data <https://example.org/a.csv>'     # a link on a searched dataset
    cols = header.split('\t')
    disp = dict(zip(['dataset_id', 'filename', 'geo_id', 'pubmed_id', 'title', 'description', 'first_author',
                     'authors', 'journal', 'year', 'tags', 'assay', 'value_type', 'species', 'display_only',
                     'links'], DISPLAY.rstrip('\n').split('\t')))
    rows.append('\t'.join(disp.get(c, '') for c in cols))
    tsv.write_text('\n'.join([header] + rows) + '\n')


@pytest.fixture(scope='module')
def client(tmp_path_factory):
    from fastapi.testclient import TestClient

    from spell_next import api
    raw = make_toy_raw(tmp_path_factory.mktemp('raw'))
    add_display_rows(raw)
    out = tmp_path_factory.mktemp('out')
    subprocess.run([sys.executable, '-m', 'spell_next.build', str(raw), str(out), '--workers', '2'],
                   check=True, capture_output=True, text=True)
    saved = api.DATA_DIR, api.LOG_DIR
    api.DATA_DIR, api.LOG_DIR = Path(out), str(tmp_path_factory.mktemp('logs'))
    with TestClient(api.app) as c:
        yield c, raw
    api.DATA_DIR, api.LOG_DIR = saved


def test_display_only_dataset_is_listed_with_links(client):
    c, _ = client
    datasets = {d['id']: d for d in c.get('/api/datasets').json()}
    assert datasets[1001]['display_only'] == 1 and datasets[1]['display_only'] == 0
    d = c.get('/api/datasets/1001').json()
    assert d['display_only'] == 1 and d['conditions'] == [] and d['species'] == 'Alpha worm'
    assert d['links'] == [{'label': 'GEO Record GSE77', 'url': 'https://example.org/geo/GSE77'},
                          {'label': 'Lab Atlas', 'url': 'https://atlas.example.org/'}]
    assert c.get('/api/datasets/1').json()['links'] == [{'label': 'Paper data', 'url': 'https://example.org/a.csv'}]
    assert c.get('/api/datasets/1001/download').status_code == 404
    info = c.get('/api/info').json()
    assert (info['datasets'], info['datasets_display_only']) == (4, 1)
    assert (info['geo_series'], info['publications']) == (4, 0)       # display-only GSE77 is not counted


def test_display_only_datasets_are_never_searched(client):
    c, _ = client
    body = c.post('/api/search', json={'genes': ['a-1', 'a-2'], 'tags': ['stress']}).json()
    assert body['datasets_searched'] == 2 and 1001 not in [d['id'] for d in body['datasets']]
    assert c.post('/api/search', json={'genes': ['a-1'], 'dataset_ids': [1001]}).status_code == 422
    assert c.post('/api/search', json={'genes': ['a-1'], 'assay': 'Single-Cell RNA-Seq'}).status_code == 422


def test_validate_accepts_display_only_rows_and_reserves_their_ids(client):
    _, raw = client
    report = validate.validate(str(raw), workers=2)
    assert not report.errors
    assert not [i for i in report.issues if 'name files not in dataset_list' in i.message]
    assert report.next_dataset_id == 1002
