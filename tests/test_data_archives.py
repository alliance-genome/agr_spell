"""Downloadable expression archives: one .tar.gz per assay, microarrays split by channel count."""
import csv
import io
import json
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

from conftest import make_toy_raw
from spell_next import validate


def add_channels(raw, channels):
    """Add a channels column to the toy dataset table: {filename: '1' or '2'}."""
    path = raw / 'meta' / 'datasets.tsv'
    rows = list(csv.DictReader(path.open(), delimiter='\t'))
    cols = list(rows[0]) + ['channels']
    with path.open('w', newline='') as fh:
        w = csv.DictWriter(fh, cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows({**r, 'channels': channels.get(r['filename'], '')} for r in rows)


@pytest.fixture(scope='module')
def built(tmp_path_factory):
    raw = make_toy_raw(tmp_path_factory.mktemp('raw'))
    (raw / 'mod.toml').write_text((raw / 'mod.toml').read_text() + '\n[build]\ndata_archives = true\n')
    add_channels(raw, {'A_set1.pcl': '2'})             # B_set1.pcl, the other microarray, has no channel count
    out = tmp_path_factory.mktemp('out')
    subprocess.run([sys.executable, '-m', 'spell_next.build', str(raw), str(out), '--workers', '2'],
                   check=True, capture_output=True, text=True)
    return raw, out


@pytest.fixture(scope='module')
def client(built, tmp_path_factory):
    from fastapi.testclient import TestClient

    from spell_next import api
    saved = api.DATA_DIR, api.LOG_DIR
    api.DATA_DIR, api.LOG_DIR = Path(built[1]), str(tmp_path_factory.mktemp('logs'))
    with TestClient(api.app) as c:
        yield c
    api.DATA_DIR, api.LOG_DIR = saved


def test_one_archive_per_assay_and_channel_count(built):
    raw, out = built
    manifest = json.loads((out / 'downloads' / 'manifest.json').read_text())
    assert [(a['file'], a['datasets']) for a in manifest] == [
        ('tm_2-channel_microarray.tar.gz', 1), ('tm_microarray.tar.gz', 1), ('tm_RNA-seq.tar.gz', 2)]
    with tarfile.open(out / 'downloads' / 'tm_RNA-seq.tar.gz') as tar:
        assert sorted(tar.getnames()) == ['tm_RNA-seq/A_set0.pcl', 'tm_RNA-seq/B_set0.pcl', 'tm_RNA-seq/datasets.tsv']
        assert tar.extractfile('tm_RNA-seq/A_set0.pcl').read() == (raw / 'expression' / 'A_set0.pcl').read_bytes()
        table = list(csv.DictReader(io.StringIO(tar.extractfile('tm_RNA-seq/datasets.tsv').read().decode()),
                                    delimiter='\t'))
    assert [(r['filename'], r['geo_id'], r['value_type'], r['conditions']) for r in table] == [
        ('A_set0.pcl', 'GSE9001', 'log2 TPM', '10'), ('B_set0.pcl', 'GSE9003', 'log2 TPM', '10')]


def test_archives_are_listed_and_served(client):
    info = client.get('/api/info').json()
    assert [a['label'] for a in info['data_archives']] == ['2-channel microarray', 'microarray', 'RNA-seq']
    assert info['built']
    r = client.get('/api/downloads/data/tm_2-channel_microarray.tar.gz')
    assert r.status_code == 200 and r.headers['content-type'] == 'application/gzip'
    with tarfile.open(fileobj=io.BytesIO(r.content)) as tar:
        assert 'tm_2-channel_microarray/A_set1.pcl' in tar.getnames()
    assert client.get('/api/downloads/data/other.tar.gz').status_code == 404
    assert client.get('/api/downloads/data/..%2Fspell.sqlite').status_code == 404


def test_validate_warns_about_microarrays_without_channels(built):
    raw, _ = built
    messages = [i.message for i in validate.validate(str(raw), 1).warnings]
    assert any(m.startswith('1 microarray datasets have no channels') and 'B_set1.pcl' in m for m in messages)


def test_no_archives_unless_configured(tmp_path):
    raw = make_toy_raw(tmp_path / 'raw')
    out = tmp_path / 'out'
    subprocess.run([sys.executable, '-m', 'spell_next.build', str(raw), str(out), '--workers', '1'],
                   check=True, capture_output=True, text=True)
    assert not (out / 'downloads').exists()
