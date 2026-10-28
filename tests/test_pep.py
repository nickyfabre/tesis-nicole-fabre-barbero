import numpy as np
import pytest

from taxoshare.pep import estimate_peps, decoy_consistency, pava_nondecreasing


def test_pava_is_monotone_and_preserves_mean():
    y = np.array([3.0, 1.0, 2.0, 5.0, 4.0])
    f = pava_nondecreasing(y)
    assert np.all(np.diff(f) >= -1e-12)
    assert f.mean() == pytest.approx(y.mean())


def test_pava_leaves_monotone_input_untouched():
    y = np.array([1.0, 2.0, 3.0])
    assert np.allclose(pava_nondecreasing(y), y)


def _competitive_tda(rng, n_spec=40000, pi0=0.25):
    n_bad = int(pi0 * n_spec)
    good = rng.normal(6.0, 1.0, n_spec - n_bad)
    bad_decoy = rng.random(n_bad) < 0.5
    bad = rng.normal(2.0, 1.0, n_bad)
    scores = np.concatenate([good, bad])
    is_decoy = np.concatenate([np.zeros(n_spec - n_bad, bool), bad_decoy])
    return scores, is_decoy


def test_peps_are_probabilities_and_monotone_in_score():
    rng = np.random.default_rng(7)
    scores, is_decoy = _competitive_tda(rng)
    peps = estimate_peps(scores, is_decoy)
    assert peps.min() >= 0.0 and peps.max() <= 1.0
    order = np.argsort(-scores)
    assert np.all(np.diff(peps[order]) >= -1e-9)


def test_pep_sum_matches_decoy_count_under_competition():
    rng = np.random.default_rng(11)
    scores, is_decoy = _competitive_tda(rng)
    peps = estimate_peps(scores, is_decoy)
    total, decoys = decoy_consistency(peps, is_decoy)
    assert total == pytest.approx(decoys, rel=0.05)


def test_ties_share_a_pep():
    scores = np.array([5.0, 5.0, 1.0, 1.0])
    is_decoy = np.array([False, True, False, True])
    peps = estimate_peps(scores, is_decoy)
    assert peps[0] == pytest.approx(peps[1])


def test_requires_decoys():
    with pytest.raises(ValueError):
        estimate_peps(np.array([1.0, 2.0]), np.array([False, False]))
