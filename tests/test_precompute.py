"""The all-genes precompute must equal one search per gene."""
import numpy as np
import pytest

from spell_next.compendium import Compendium, Dataset
from spell_next.precompute import precompute_all
from spell_next.search import search


def _compendium(seed=3):
    rng = np.random.default_rng(seed)
    n_genes = 30
    datasets = []
    for d in range(8):
        n_conds = 4 if d == 7 else 9                       # one dataset too short to search
        base = rng.normal(size=n_conds)
        values = rng.normal(size=(n_genes, n_conds))
        values[:6] = base + 0.2 * rng.normal(size=(6, n_conds))      # a co-expressed module
        values[6] = 5.0 + 1e-9 * rng.normal(size=n_conds)            # near-constant gene
        values[7] = np.where(np.arange(n_conds) == 2, 0.4, 0.0)      # identical "spike" genes
        values[8] = np.where(np.arange(n_conds) == 2, 0.2, 0.0)
        values[rng.random(values.shape) < 0.08] = np.nan             # missing values
        keep = np.sort(rng.choice(n_genes, n_genes - d, replace=False))   # genes vary by dataset
        datasets.append(Dataset(100 + d, 'd%d.pcl' % d, ['c%d' % i for i in range(n_conds)],
                                keep.astype(np.int32), values[keep], z_mean=0.05, z_std=0.5))
    return Compendium(list(range(1000, 1000 + n_genes)), ['G%d' % i for i in range(n_genes)], datasets)


def test_precompute_matches_per_gene_search(tmp_path):
    comp = _compendium()
    npz = str(tmp_path / 'c.npz')
    comp.save(npz)
    comp = Compendium.load(npz)
    results = precompute_all(comp, npz, top_genes=10, top_datasets=10, workers=2, log=lambda *_: None)
    assert set(results) == set(comp.gene_ids)
    for gid in comp.gene_ids:
        r = search(comp, [gid], max_genes=10, precision='float64')
        hits = [(g, s) for g, s in zip(r.gene_ids, r.gene_scores) if g != gid][:10]
        tot = sum(r.raw_weights) or 1.0
        dsets = [(d, w / tot * 100) for d, w in zip(r.dataset_ids, r.raw_weights)][:10]
        got_hits, got_dsets = results[gid]
        assert [g for g, _ in got_hits] == [g for g, _ in hits], gid
        assert [s for _, s in got_hits] == pytest.approx([s for _, s in hits], abs=1e-9)
        assert [d for d, _ in got_dsets] == [d for d, _ in dsets], gid
        assert [w for _, w in got_dsets] == pytest.approx([w for _, w in dsets], abs=1e-9)
