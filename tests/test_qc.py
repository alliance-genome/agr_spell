import numpy as np
import pytest

from spell_next import qc
from spell_next.search import _Arith, pearson_with_row


def _diag(values):
    return qc.diagnose(np.asarray(values, dtype=float))


def test_pairwise_pearson_matches_row_by_row_with_missing_values():
    rng = np.random.default_rng(3)
    values = rng.normal(size=(15, 9))
    values[rng.random(values.shape) < 0.15] = np.nan
    full = qc.pairwise_pearson(values)
    for i in (0, 4, 11):
        np.testing.assert_allclose(full[i], pearson_with_row(values, i, _Arith('float64')), atol=1e-12)


def test_pair_stats_of_identical_rows_is_large():
    base = np.arange(10, dtype=float)
    z_mean, z_std = qc.pair_stats(np.vstack([base, base * 2, base + 3]))
    assert z_mean > 7 and z_std == pytest.approx(0, abs=1e-6)


@pytest.mark.parametrize('values, expected', [
    ([[-1.2, 0.3, 2.1], [0.4, -0.8, 1.5]], 'log'),                       # log ratios
    ([[8.1, 9.3, 12.0], [7.7, 10.2, 11.4]], 'log'),                       # log intensities
    ([[0.3, 1.0, 2.5], [0.7, 1.4, 6.2], [0.5, 0.9, 1.1]], 'ratio'),       # unlogged ratios
    ([[120.0, 3500.0, 40.0], [800.0, 15.0, 9000.0]], 'intensity'),        # raw counts
    ([[1.2, 1.8, 5.0], [1.5, 2.2, 1.9]], 'unclear'),                      # positive, never below 1
])
def test_classify_scale(values, expected):
    assert qc.classify_scale(_diag(values)) == expected


def test_normalize_logs_ratios_and_does_not_centre_by_default():
    values = np.array([[0.5, 1.0, 2.0], [0.25, 1.0, 8.0], [1.0, 1.0, 4.0]])
    out, actions, warnings = qc.normalize(values)
    np.testing.assert_allclose(out[0], [-1.0, 0.0, 1.0])
    assert any('log2 transform' in a for a in actions)
    assert not any('centred' in a for a in actions)


def test_normalize_can_centre_conditions():
    values = np.array([[1.0, 5.0], [2.0, 6.0], [3.0, 7.0]])
    out, actions, _ = qc.normalize(values, centre=True)
    np.testing.assert_allclose(np.median(out, axis=0), [0.0, 0.0])
    assert any('centred' in a for a in actions)


def test_normalize_flags_shifted_conditions():
    values = np.array([[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]])
    _, _, warnings = qc.normalize(values)
    assert any('shifted' in w for w in warnings)
