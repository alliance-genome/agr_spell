"""End-to-end build of a small two-species MOD from generic gene and dataset tables."""
import gzip
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from spell_next import modconfig

SPECIES = {'A': 'Alpha worm', 'B': 'Beta worm'}
GENES_PER_SPECIES = 12
N_CONDS = 10

CONFIG = '''
[mod]
id = "TM"
name = "Test MOD"
organism = "test"
examples = ["a-1", "a-1 a-2"]

[links]
gene = "https://example.org/gene/{mod_id}"

[genes]
file = "meta/genes.tsv"
format = "gene_table"

[datasets]
metadata = "meta/datasets.tsv"
format = "dataset_table"

[go]
obo = "meta/go.obo"
gaf = ["meta/a.gaf.gz", "meta/b.gaf"]
'''

OBO = '''format-version: 1.2
data-version: releases/test

[Term]
id: GO:0000001
name: everything
namespace: biological_process

[Term]
id: GO:0000002
name: alpha module
namespace: biological_process
is_a: GO:0000001 ! everything

[Term]
id: GO:0000003
name: beta module
namespace: biological_process
is_a: GO:0000001 ! everything
'''


def gaf_line(mod_id, term):
    return '\t'.join(['TM', mod_id, mod_id, '', term, 'PMID:1', 'IDA', '', 'P', '', '', 'gene', 'taxon:1',
                      '20260101', 'TM', '', '']) + '\n'


@pytest.fixture(scope='module')
def build(tmp_path_factory):
    raw = tmp_path_factory.mktemp('raw')
    out = tmp_path_factory.mktemp('out')
    (raw / 'expression').mkdir()
    (raw / 'meta').mkdir()
    (raw / 'mod.toml').write_text(CONFIG)
    (raw / 'meta' / 'go.obo').write_text(OBO)

    gene_list, gene_rows = [], ['systematic\tmod_id\tstandard\taliases\tspecies\tdescription']
    gid = 1000
    for sp, species in SPECIES.items():
        for i in range(1, GENES_PER_SPECIES + 1):
            gene_list.append('%d\t%s%d' % (gid, sp, i))
            gene_rows.append('%s%d\tMOD%s%d\t%s-%d\t%s%d-alias\t%s\tgene %d of %s'
                             % (sp, i, sp, i, sp.lower(), i, sp.lower(), i, species, i, species))
            gid += 1
    (raw / 'gene_list.txt').write_text('\n'.join(gene_list) + '\n')
    (raw / 'meta' / 'genes.tsv').write_text('\n'.join(gene_rows) + '\n')

    # Two datasets per species; genes 1-4 of each species rise and fall together.
    rng = np.random.default_rng(7)
    ds_rows, ds_list = ['filename\tgeo_id\ttitle\tfirst_author\tyear\ttags'], []
    ds_id = 1
    for sp in SPECIES:
        for k in range(2):
            name = '%s_set%d.pcl' % (sp, k)
            base = rng.normal(size=N_CONDS)
            lines = ['ID\tNAME\tGWEIGHT\t' + '\t'.join('c%d' % c for c in range(N_CONDS)),
                     'EWEIGHT\t\t\t' + '\t'.join(['1'] * N_CONDS)]
            for i in range(1, GENES_PER_SPECIES + 1):
                vals = (base + 0.1 * rng.normal(size=N_CONDS)) if i <= 4 else rng.normal(size=N_CONDS)
                lines.append('%s%d\t%s%d\t1\t' % (sp, i, sp, i) + '\t'.join('%.4f' % (8 + v) for v in vals))
            (raw / 'expression' / name).write_text('\n'.join(lines) + '\n')
            (raw / 'expression' / (name + '.sinfo')).write_text('All\t0\t1\n')
            ds_list.append('%d\t%s' % (ds_id, name))
            ds_rows.append('%s\tGSE%d\t%s study %d\tAuthor %s\t2020\t%s' % (
                name, 9000 + ds_id, SPECIES[sp], k, sp, 'stress' if k == 0 else 'growth'))
            ds_id += 1
    (raw / 'dataset_list.txt').write_text('\n'.join(ds_list) + '\n')
    (raw / 'meta' / 'datasets.tsv').write_text('\n'.join(ds_rows) + '\n')

    # Species A's co-expressed genes share "alpha module"; species B's share "beta module".
    with gzip.open(raw / 'meta' / 'a.gaf.gz', 'wt') as fh:
        fh.write('!date-generated: 20260101\n')
        fh.writelines(gaf_line('MODA%d' % i, 'GO:0000002' if i <= 4 else 'GO:0000001')
                      for i in range(1, GENES_PER_SPECIES + 1))
    (raw / 'meta' / 'b.gaf').write_text('!date-generated: 20260102\n' + ''.join(
        gaf_line('MODB%d' % i, 'GO:0000003' if i <= 4 else 'GO:0000001') for i in range(1, GENES_PER_SPECIES + 1)))

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


def test_bad_config_is_rejected(tmp_path):
    (tmp_path / 'mod.toml').write_text('[genes]\nformat = "excel"\n')
    with pytest.raises(modconfig.ConfigError):
        modconfig.load(str(tmp_path))
    (tmp_path / 'mod.toml').write_text('[links]\ngene = "https://x.org/{sgdid}"\n')
    with pytest.raises(modconfig.ConfigError):
        modconfig.load(str(tmp_path))
