import numpy as np
import pytest

from taxoshare.allocate import allocate_species, nodes_from_species, rigid_lca_nodes
from taxoshare.tree import flat_two_level_tree

S = ["A", "B", "C"]


def fs(*names):
    return frozenset(names)


def test_pi_is_a_distribution():
    r = allocate_species([fs("A"), fs("B"), fs("A", "B")], np.array([10.0, 30.0, 60.0]), S)
    assert r.pi.min() >= 0
    assert r.pi.sum() == pytest.approx(1.0)


def test_em_loglik_is_monotone():
    rng = np.random.default_rng(3)
    sets, ints = [], []
    for _ in range(400):
        k = rng.integers(1, 4)
        sets.append(frozenset(rng.choice(S, size=k, replace=False).tolist()))
        ints.append(float(rng.lognormal(0, 1)))
    r = allocate_species(sets, np.array(ints), S, method="em", max_iter=500)
    trace = np.array(r.loglik_trace)
    assert np.all(np.diff(trace) >= -1e-9), "EM must not decrease the log-likelihood"
    assert r.converged


def test_two_species_shared_mass_follows_unique_evidence():
    # A unique 10, B unique 30, 60 shared: the fixed point is the unique ratio.
    r = allocate_species([fs("A"), fs("B"), fs("A", "B")], np.array([10.0, 30.0, 60.0]), ["A", "B"],
                         method="em")
    assert r.pi[0] == pytest.approx(0.25, abs=1e-6)
    assert r.pi[1] == pytest.approx(0.75, abs=1e-6)


def test_unique_arm_ignores_shared_peptides():
    r = allocate_species([fs("A"), fs("B"), fs("A", "B")], np.array([10.0, 30.0, 900.0]), ["A", "B"],
                         method="unique")
    assert r.pi[0] == pytest.approx(0.25)


def test_proportional_is_one_em_step_from_unique_init():
    sets = [fs("A"), fs("B"), fs("C"), fs("A", "B"), fs("B", "C"), fs("A", "C")]
    ints = np.array([10.0, 20.0, 30.0, 100.0, 100.0, 100.0])
    prop = allocate_species(sets, ints, S, method="proportional")
    em1 = allocate_species(sets, ints, S, method="em", init="unique", max_iter=1)
    assert np.allclose(prop.pi, em1.pi)


def test_em_improves_on_one_step_when_they_differ():
    sets = [fs("A"), fs("B"), fs("C"), fs("A", "B"), fs("B", "C"), fs("A", "C")]
    ints = np.array([10.0, 20.0, 30.0, 100.0, 100.0, 100.0])
    prop = allocate_species(sets, ints, S, method="proportional")
    em = allocate_species(sets, ints, S, method="em", init="unique")
    assert not np.allclose(prop.pi, em.pi), "this configuration should separate them"
    assert em.loglik > prop.loglik


def test_pep_of_one_removes_an_observation():
    sets = [fs("A"), fs("B"), fs("A", "B"), fs("A")]
    ints = np.array([10.0, 30.0, 60.0, 500.0])
    peps = np.array([0.0, 0.0, 0.0, 1.0])  # the 500 is a certain false identification
    with_pep = allocate_species(sets, ints, ["A", "B"], method="em", peps=peps)
    without = allocate_species(sets[:3], ints[:3], ["A", "B"], method="em")
    assert np.allclose(with_pep.pi, without.pi, atol=1e-9)


def test_pep_downweights_partially():
    sets = [fs("A"), fs("B"), fs("A")]
    ints = np.array([10.0, 10.0, 10.0])
    half = allocate_species(sets, ints, ["A", "B"], method="em", peps=np.array([0.0, 0.0, 0.5]))
    assert half.pi[0] == pytest.approx(15.0 / 25.0)


def test_uniform_init_can_reach_a_species_with_no_unique_peptides():
    # C has no unique peptide; a single proportional step can never give it signal.
    sets = [fs("A"), fs("B"), fs("B", "C"), fs("C", "A")]
    ints = np.array([10.0, 10.0, 500.0, 500.0])
    prop = allocate_species(sets, ints, S, method="proportional")
    em = allocate_species(sets, ints, S, method="em", init="uniform")
    assert prop.pi[S.index("C")] < 1e-9   # only the initialisation floor, no real signal
    assert em.pi[S.index("C")] > 0.05


def test_nodes_roll_up_and_rigid_lca_is_normalised():
    tree = flat_two_level_tree({"g1": ["A", "B"], "g2": ["C"]})
    node_vals = nodes_from_species(tree, {"A": 0.2, "B": 0.3, "C": 0.5})
    assert node_vals["g1"] == pytest.approx(0.5)
    assert node_vals["root"] == pytest.approx(1.0)

    rigid = rigid_lca_nodes(tree, ["A", "g1", "C"], np.array([10.0, 40.0, 50.0]))
    assert rigid["root"] == pytest.approx(1.0)
    assert rigid["g1"] == pytest.approx(0.5)   # 10 at A + 40 parked at g1
    assert rigid["A"] == pytest.approx(0.1)


def test_all_rank_candidate_sets_span_lca_down_to_members():
    from taxoshare.allocate import all_rank_candidate_sets, node_totals_at_or_below
    tree = flat_two_level_tree({"g1": ["A", "B"], "g2": ["C"]})
    sets = all_rank_candidate_sets(tree, [fs("A"), fs("A", "B"), fs("A", "C")])
    assert sets[0] == fs("A")
    assert sets[1] == fs("A", "B", "g1")
    assert sets[2] == fs("A", "C", "g1", "g2", "root")


def test_node_totals_at_or_below_accumulate_upwards():
    from taxoshare.allocate import node_totals_at_or_below
    tree = flat_two_level_tree({"g1": ["A", "B"], "g2": ["C"]})
    tot = node_totals_at_or_below(tree, {"A": 0.3, "g1": 0.2, "C": 0.5})
    assert tot["g1"] == pytest.approx(0.5)   # 0.3 at A plus 0.2 held at the genus
    assert tot["root"] == pytest.approx(1.0)
