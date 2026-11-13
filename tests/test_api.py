"""API tests against a real build in data/build (skipped when it is absent)."""
import os
from pathlib import Path

import pytest

BUILD = Path(__file__).resolve().parent.parent / 'data' / 'build'
pytestmark = pytest.mark.skipif(not (BUILD / 'spell.sqlite').exists(), reason='no data/build')


@pytest.fixture(scope='module')
def log_dir(tmp_path_factory):
    return tmp_path_factory.mktemp('logs')


@pytest.fixture(scope='module')
def client(log_dir):
    os.environ['SPELL_NEXT_DATA'] = str(BUILD)
    from fastapi.testclient import TestClient

    from spell_next import api
    api.DATA_DIR = BUILD
    api.LOG_DIR = str(log_dir)       # keep test searches out of the real logs
    with TestClient(api.app) as c:
        yield c


def test_lookup_resolves_standard_and_systematic_names(client):
    body = client.get('/api/genes/lookup', params={'q': 'act1, YAL001C nosuchgene'}).json()
    assert [g['systematic'] for g in body['genes']] == ['YFL039C', 'YAL001C']
    assert body['unknown'] == ['nosuchgene']


def test_precomputed_single_gene_result(client):
    act1 = client.get('/api/genes/lookup', params={'q': 'ACT1'}).json()['genes'][0]
    body = client.get('/api/genes/%d/coexpression' % act1['id'], params={'limit': 20}).json()
    assert body['source'] == 'precomputed'
    assert len(body['genes']) == 20
    assert all(g['id'] != act1['id'] for g in body['genes'])
    scores = [g['score'] for g in body['genes']]
    assert scores == sorted(scores, reverse=True)
    assert body['datasets'][0]['weight'] >= body['datasets'][-1]['weight']


def test_live_multi_gene_search(client):
    body = client.post('/api/search', json={'genes': ['RPL3', 'RPL4A', 'RPS2'], 'limit': 10}).json()
    assert body['source'] == 'live'
    assert len(body['query']) == 3
    assert len(body['genes']) == 10
    assert sum(1 for g in body['genes'] if (g['standard'] or '').startswith('RP')) >= 7
    assert sum(d['weight'] for d in body['datasets']) <= 100.0 + 1e-6


def test_search_with_no_known_genes_is_rejected(client):
    assert client.post('/api/search', json={'genes': ['nosuchgene']}).status_code == 422


def test_dataset_detail_and_expression(client):
    act1 = client.get('/api/genes/lookup', params={'q': 'ACT1'}).json()['genes'][0]
    ds_id = client.get('/api/genes/%d/coexpression' % act1['id']).json()['datasets'][0]['id']
    detail = client.get('/api/datasets/%d' % ds_id).json()
    assert len(detail['conditions']) == detail['num_conds']
    expr = client.get('/api/datasets/%d/expression' % ds_id, params={'genes': str(act1['id'])}).json()
    assert len(expr['rows'][0]['values']) == detail['num_conds']


def test_front_end_is_served(client):
    assert 'SPELL' in client.get('/').text
    assert client.get('/static/app.js').status_code == 200


def test_go_enrichment_of_ribosomal_genes(client):
    names = 'RPL3 RPL4A RPL5 RPL6A RPL7A RPS2 RPS3 RPS5 RPS6A RPS8A'
    ids = [g['id'] for g in client.get('/api/genes/lookup', params={'q': names}).json()['genes']]
    body = client.post('/api/go-enrichment', json={'gene_ids': ids, 'aspect': 'P'}).json()
    top = body['terms'][0]
    assert 'translation' in top['name'] or 'ribosom' in top['name']
    assert top['list_count'] == 10 and top['fdr'] < 1e-6
    assert body['versions']['go_ontology']


def test_expression_profile_with_raw_values(client):
    act1 = client.get('/api/genes/lookup', params={'q': 'ACT1'}).json()['genes'][0]
    ds_id = client.get('/api/genes/%d/coexpression' % act1['id']).json()['datasets'][0]['id']
    row = client.get('/api/datasets/%d/expression' % ds_id,
                     params={'genes': str(act1['id']), 'raw': 'true'}).json()['rows'][0]
    present = [(c, r) for c, r in zip(row['values'], row['raw']) if c is not None]
    offsets = {round(r - c, 2) for c, r in present}
    assert len(offsets) == 1                     # centred = raw - the gene's mean
    assert abs(sum(c for c, _ in present)) < 1e-2 * len(present)


def test_dataset_listing_and_tags(client):
    datasets = client.get('/api/datasets').json()
    tags = client.get('/api/tags').json()
    assert len(datasets) > 700
    stress = next(t for t in tags if t['tag'] == 'stress')
    assert stress['count'] == sum('stress' in d['tags'] for d in datasets)


def test_search_restricted_to_tags(client):
    tagged = {d['id'] for d in client.get('/api/datasets').json() if 'heat shock' in d['tags']}
    body = client.post('/api/search', json={'genes': ['HSP104'], 'tags': ['heat shock']}).json()
    assert body['source'] == 'live' and body['datasets_searched'] == len(tagged)
    assert {d['id'] for d in body['datasets']} <= tagged
    assert client.post('/api/search', json={'genes': ['HSP104'], 'tags': ['no such tag']}).status_code == 422


def test_expression_listing(client):
    body = client.get('/api/expression-listing', params={'genes': 'GAL1 GAL7 GAL10'}).json()
    assert [g['standard'] for g in body['query']] == ['GAL1', 'GAL7', 'GAL10']
    measured = [d['genes_measured'] for d in body['datasets']]
    assert measured == sorted(measured, reverse=True) and measured[0] == 3
    tagged = client.get('/api/expression-listing', params={'genes': 'GAL1', 'tags': 'heat shock'}).json()
    assert 0 < len(tagged['datasets']) < len(body['datasets'])


def test_dataset_download(client):
    act1 = client.get('/api/genes/lookup', params={'q': 'ACT1'}).json()['genes'][0]
    ds_id = client.get('/api/genes/%d/coexpression' % act1['id']).json()['datasets'][0]['id']
    detail = client.get('/api/datasets/%d' % ds_id).json()
    res = client.get('/api/datasets/%d/download' % ds_id)
    assert res.status_code == 200 and 'attachment' in res.headers['content-disposition']
    lines = res.text.rstrip('\n').split('\n')
    assert lines[0].split('\t')[:2] == ['systematic_name', 'standard_name']
    assert len(lines[0].split('\t')) == detail['num_conds'] + 2
    assert len(lines) == detail['num_genes'] + 1
    assert any(line.startswith('YFL039C\tACT1\t') for line in lines)


def test_embed_page_and_info(client):
    assert 'embed.js' in client.get('/embed', params={'gene': 'ACT1'}).text
    info = client.get('/api/info').json()
    assert info['datasets'] > 700 and info['versions']['go_ontology']


def test_searches_are_logged_without_personal_data(client, log_dir):
    client.post('/api/search', json={'genes': ['CLN1', 'CLN2'], 'limit': 5})
    lines = [ln for f in log_dir.glob('searches_*.tsv') for ln in f.read_text().splitlines()]
    multi = [ln.split('\t') for ln in lines if '\tmulti\t' in ln]
    assert multi and multi[-1][3] == 'YMR199W,YPL256C'
    assert all(len(ln.split('\t')) == 5 for ln in lines)
