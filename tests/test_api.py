"""API tests against a real build in data/build (skipped when it is absent)."""
import os
from pathlib import Path

import pytest

BUILD = Path(__file__).resolve().parent.parent / 'data' / 'build'
pytestmark = pytest.mark.skipif(not (BUILD / 'spell.sqlite').exists(), reason='no data/build')


@pytest.fixture(scope='module')
def client():
    os.environ['SPELL_NEXT_DATA'] = str(BUILD)
    from fastapi.testclient import TestClient

    from spell_next import api
    api.DATA_DIR = BUILD
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
