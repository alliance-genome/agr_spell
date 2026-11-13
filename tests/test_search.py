import math

import numpy as np
import pytest

from spell_next.compendium import Compendium, Dataset, parse_float, read_pcl, tokenize
from spell_next.search import _Arith, fisher_z, pearson_with_row, search, zscore2


def test_tokenize_collapses_tabs_like_spell():
    assert tokenize('A\t\t1\t2', faithful=True) == ['A', '1', '2']
    assert tokenize('A\t\t1\t2', faithful=False) == ['A', '', '1', '2']


@pytest.mark.parametrize('token, expected', [
    ('1.5', 1.5), (' -2e3', -2000.0), ('3.2abc', 3.2), ('.5', 0.5),
])
def test_parse_float_numeric_prefix(token, expected):
    assert parse_float(token) == expected


@pytest.mark.parametrize('token', ['', 'NA', 'nan', 'abc'])
def test_parse_float_missing(token):
    assert math.isnan(parse_float(token))


def test_read_pcl_with_gweight_and_eweight(tmp_path):
    pcl = tmp_path / 'x.pcl'
    pcl.write_text('YORF\tNAME\tGWEIGHT\tc1\tc2\tc3\n'
                   'EWEIGHT\t\t\t1\t1\t1\n'
                   'YAL001C\tTFC3\t1\t1.0\t2.0\t3.0\n'
                   'unknown\tACT1\t1\t4.0\tNA\t6.0\n'
                   'NOTAGENE\tNOPE\t1\t1\t1\t1\n')
    conds, gene_index, values = read_pcl(str(pcl), {'YAL001C': 0, 'ACT1': 1})
    assert conds == ['c1', 'c2', 'c3']
    assert list(gene_index) == [0, 1]
    np.testing.assert_allclose(values[0], [1.0, 2.0, 3.0])
    assert math.isnan(values[1, 1])


def test_pearson_matches_numpy_on_complete_pairs():
    rng = np.random.default_rng(0)
    values = rng.normal(size=(20, 12))
    values[3, 4] = np.nan
    rho = pearson_with_row(values, 0, _Arith('float64'))
    ok = ~np.isnan(values[3])
    assert rho[0] == pytest.approx(1.0)
    assert rho[3] == pytest.approx(np.corrcoef(values[0, ok], values[3, ok])[0, 1])


def test_fisher_z_clips_perfect_correlation():
    z = fisher_z(np.array([1.0, -1.0, 0.0]), _Arith('float64'))
    assert z[0] == pytest.approx(0.5 * math.log(1.999999 / 0.000001))
    assert z[1] == pytest.approx(-z[0])
    assert z[2] == 0.0


def test_zscore2_zeroes_below_one_and_squares_above():
    ds = Dataset(1, 'x', [], np.array([], dtype=np.int32), np.zeros((0, 0)), z_mean=0.0, z_std=1.0)
    rho = np.tanh(np.array([0.5, 2.0, np.nan]))      # Fisher z of tanh(x) is x
    out = zscore2(ds, rho, _Arith('float64'))
    assert out[0] == 0.0
    assert out[1] == pytest.approx(4.0)
    assert math.isnan(out[2])


def _toy_compendium(num_datasets=6, num_genes=8, num_conds=10, seed=1):
    """Gene 1 tracks gene 0 in every dataset; the other genes are noise."""
    rng = np.random.default_rng(seed)
    datasets = []
    for d in range(num_datasets):
        base = rng.normal(size=num_conds)
        values = rng.normal(size=(num_genes, num_conds))
        values[0] = base
        values[1] = base + 0.05 * rng.normal(size=num_conds)
        datasets.append(Dataset(100 + d, 'd%d.pcl' % d, ['c%d' % i for i in range(num_conds)],
                                np.arange(num_genes, dtype=np.int32), values, z_mean=0.0, z_std=0.5))
    return Compendium(list(range(1000, 1000 + num_genes)), ['G%d' % i for i in range(num_genes)], datasets)


@pytest.mark.parametrize('precision', ['float32', 'float64'])
def test_single_gene_query_ranks_itself_then_partner(precision):
    result = search(_toy_compendium(), [1000], precision=precision, max_genes=3)
    assert result.gene_ids[:2] == [1000, 1001]
    assert result.dataset_weights is None
    assert len(result.gene_ids) == 4                  # max_genes plus the query gene
    assert sorted(result.dataset_ids) == list(range(100, 106))


def test_multi_gene_query_weights_sum_to_100():
    result = search(_toy_compendium(), [1000, 1001], precision='float64')
    assert set(result.gene_ids[:2]) == {1000, 1001}
    assert sum(result.dataset_weights) == pytest.approx(100.0)


def test_datasets_with_fewer_than_five_conditions_are_skipped():
    comp = _toy_compendium(num_conds=4)
    result = search(comp, [1000, 1001], precision='float64')
    assert result.dataset_ids == []
