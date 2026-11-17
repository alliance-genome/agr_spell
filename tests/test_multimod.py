"""End-to-end build of a small two-species MOD from generic gene and dataset tables."""
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import GENES_PER_SPECIES, make_toy_raw
from spell_next import modconfig

@pytest.fixture(scope='module')
def build(tmp_path_factory):
    raw = make_toy_raw(tmp_path_factory.mktemp('raw'))
    out = tmp_path_factory.mktemp('out')
    subprocess.run([sys.executable, '-m', 'spell_next.build', str(raw), str(out), '--workers', '2'],
                   check=True, capture_output=True, text=True)
    return out


@pytest.fixture(scope='module')
def client(build, tmp_path_factory):
    from fastapi.testclient import TestClient

    from spell_next import api
    saved = api.DATA_DIR, api.LOG_DIR
    api.DATA_DIR, api.LOG_DIR = Path(build), str(tmp_path_factory.mktemp('logs'))
    with TestClient(api.app) as c:
        yield c
    api.DATA_DIR, api.LOG_DIR = saved


def test_species_assigned_to_genes_and_datasets(build):
    db = sqlite3.connect(os.path.join(build, 'spell.sqlite'))
    assert dict(db.execute('SELECT systematic, species FROM genes WHERE systematic IN ("A1", "B1")')) == \
        {'A1': 'Alpha worm', 'B1': 'Beta worm'}
    assert dict(db.execute('SELECT filename, species FROM datasets')) == {
        'A_set0.pcl': 'Alpha worm', 'A_set1.pcl': 'Alpha worm', 'B_set0.pcl': 'Beta worm', 'B_set1.pcl': 'Beta worm'}
    assert db.execute('SELECT mod_id, standard FROM genes WHERE systematic = "B2"').fetchone() == ('MODB2', 'b-2')
    info = dict(db.execute('SELECT key, value FROM build_info'))
    assert info['mod'] == 'TM' and info['go_annotations'] == '20260101, 20260102'


def test_config_reaches_the_front_end(client):
    cfg = client.get('/api/config').json()
    assert cfg['mod']['id'] == 'TM' and cfg['mod']['examples'][0] == 'a-1'
    assert cfg['links']['gene'] == 'https://example.org/gene/{mod_id}'
    assert sorted(cfg['species']) == ['Alpha worm', 'Beta worm']


def test_hits_stay_within_the_query_species(client):
    a1 = client.get('/api/genes/lookup', params={'q': 'a-1'}).json()['genes'][0]
    body = client.get('/api/genes/%d/coexpression' % a1['id'], params={'limit': 10}).json()
    assert {g['species'] for g in body['genes'] if g['score'] > 0} == {'Alpha worm'}
    assert {g['standard'] for g in body['genes'][:3]} == {'a-2', 'a-3', 'a-4'}
    assert {d['species'] for d in body['datasets']} == {'Alpha worm'}


def test_go_enrichment_uses_the_species_background(client):
    lookup = client.get('/api/genes/lookup', params={'q': 'b-1 b-2 b-3 b-4'}).json()
    body = client.post('/api/go-enrichment', json={'gene_ids': [g['id'] for g in lookup['genes']]}).json()
    assert body['species'] == 'Beta worm'
    assert body['background_size'] == GENES_PER_SPECIES          # only species B's annotated genes
    assert body['terms'][0]['name'] == 'beta module' and body['terms'][0]['list_count'] == 4


def test_aliases_and_dataset_table_metadata(client):
    assert client.get('/api/genes/lookup', params={'q': 'b3-alias'}).json()['genes'][0]['systematic'] == 'B3'
    datasets = {d['geo_id']: d for d in client.get('/api/datasets').json()}
    assert datasets['GSE9003']['title'] == 'Beta worm study 0' and datasets['GSE9003']['tags'] == ['stress']


def test_assay_fields_and_counts(client):
    assert client.get('/api/assays').json() == [{'assay': 'RNA-seq', 'count': 2}, {'assay': 'microarray', 'count': 2}]
    datasets = {d['geo_id']: d for d in client.get('/api/datasets').json()}
    assert (datasets['GSE9001']['assay'], datasets['GSE9001']['value_type']) == ('RNA-seq', 'log2 TPM')
    assert (datasets['GSE9002']['assay'], datasets['GSE9002']['value_type']) == ('microarray', None)
    assert client.get('/api/datasets/1').json()['assay'] == 'RNA-seq'
    assert client.get('/api/config').json()['assays'] == ['microarray', 'RNA-seq']


@pytest.mark.parametrize('genes', [['a-1'], ['a-1', 'a-2']])
def test_search_restricted_to_an_assay(client, genes):
    body = client.post('/api/search', json={'genes': genes, 'assay': 'RNA-seq'}).json()
    assert body['assay'] == 'RNA-seq' and body['datasets_searched'] == 2
    assert {d['assay'] for d in body['datasets']} == {'RNA-seq'}
    both = client.post('/api/search', json={'genes': genes, 'assay': 'microarray', 'tags': ['growth']}).json()
    assert {(d['assay'], d['tags']) for d in both['datasets']} == {('microarray', 'growth')}
    assert client.post('/api/search', json={'genes': genes, 'assay': 'microarray', 'tags': ['stress']}
                       ).status_code == 422                        # no microarray dataset is tagged stress


def test_search_restricted_to_listed_datasets(client):
    ids = {d['geo_id']: d['id'] for d in client.get('/api/datasets').json()}
    body = client.post('/api/search', json={'genes': ['a-1', 'a-2'], 'dataset_ids': [ids['GSE9002']]}).json()
    assert body['datasets_searched'] == 1 and [d['geo_id'] for d in body['datasets']] == ['GSE9002']
    both = client.post('/api/search', json={'genes': ['a-1'], 'assay': 'RNA-seq',
                                            'dataset_ids': [ids['GSE9001'], ids['GSE9002']]}).json()
    assert [d['geo_id'] for d in both['datasets']] == ['GSE9001']      # ids and assay combine
    assert client.post('/api/search', json={'genes': ['a-1'], 'dataset_ids': [999]}).status_code == 422


def test_expression_listing_restricted_to_an_assay(client):
    body = client.get('/api/expression-listing', params={'genes': 'a-1', 'assay': 'microarray'}).json()
    assert body['assay'] == 'microarray' and [d['geo_id'] for d in body['datasets']] == ['GSE9002']


def test_tags_carry_keyword_categories_definitions_and_links(client):
    tags = {t['tag']: t for t in client.get('/api/tags').json()}
    assert set(tags) == {'stress', 'growth'}                      # only tags in use; heat shock is not
    assert (tags['stress']['category'], tags['stress']['definition']) == \
        ('Stress & stimuli', 'Responses to stress, of any kind.')
    assert tags['growth']['link'] == 'https://example.org/keyword/growth_kw'
    assert tags['growth']['position'] < tags['stress']['position']   # keyword-list order


def test_bad_config_is_rejected(tmp_path):
    (tmp_path / 'mod.toml').write_text('[genes]\nformat = "excel"\n')
    with pytest.raises(modconfig.ConfigError):
        modconfig.load(str(tmp_path))
    (tmp_path / 'mod.toml').write_text('[links]\ngene = "https://x.org/{sgdid}"\n')
    with pytest.raises(modconfig.ConfigError):
        modconfig.load(str(tmp_path))
