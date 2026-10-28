"""A ground-truth community for testing allocation strategies.

The generator mirrors the features that make real metaproteomic data hard:
peptide sharing that decays with phylogenetic distance, a wide dynamic range of
species abundance, species that are genuinely absent, and a population of false
identifications whose scores overlap the decoys.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .tree import Tree, flat_two_level_tree


@dataclass
class SimulatedRun:
    tree: Tree
    species: list[str]
    true_abundance: np.ndarray          # normalised, over `species`
    candidate_sets: list[frozenset[str]]  # true containing set per target PSM
    lca_nodes: list[str]                  # LCA of the containing set
    lca_descendant_sets: list[frozenset[str]]  # leaves under the LCA (what node-level sees)
    intensities: np.ndarray
    is_false: np.ndarray                # target PSM is a false identification
    scores: np.ndarray                  # all PSMs: targets then decoys
    is_decoy: np.ndarray


def make_tree(n_genera: int = 6, per_genus: int = 5) -> tuple[Tree, list[str]]:
    groups = {
        f"g{gi}": [f"g{gi}s{si}" for si in range(per_genus)]
        for gi in range(n_genera)
    }
    tree = flat_two_level_tree(groups)
    return tree, tree.leaves


def simulate(
    rng: np.random.Generator,
    *,
    n_genera: int = 6,
    per_genus: int = 5,
    n_peptides: int = 6000,
    frac_absent: float = 0.2,
    within_genus_share: float = 0.025,
    across_genus_share: float = 0.001,
    false_rate: float = 0.05,
    dynamic_range: float = 2.0,
) -> SimulatedRun:
    tree, species = make_tree(n_genera, per_genus)
    n = len(species)
    genus_of = {s: tree.parent[s] for s in species}

    abundance = rng.lognormal(0.0, dynamic_range, n)
    n_absent = int(frac_absent * n)
    if n_absent:
        abundance[rng.choice(n, size=n_absent, replace=False)] = 0.0
    if abundance.sum() <= 0:
        abundance[0] = 1.0
    true_abundance = abundance / abundance.sum()
    abund = dict(zip(species, true_abundance))

    cand, lca_nodes, lca_desc, intensities, is_false = [], [], [], [], []

    for _ in range(n_peptides):
        origin = species[rng.integers(n)]
        contains = {origin}
        for other in species:
            if other == origin:
                continue
            p = within_genus_share if genus_of[other] == genus_of[origin] else across_genus_share
            if rng.random() < p:
                contains.add(other)
        signal = sum(abund[s] for s in contains)
        if signal <= 0:
            continue  # every carrier is absent: peptide is not observed
        intensity = signal * rng.lognormal(0.0, 0.4)
        node = tree.lca(contains)
        cand.append(frozenset(contains))
        lca_nodes.append(node)
        lca_desc.append(tree.descendant_leaves(node))
        intensities.append(intensity)
        is_false.append(False)

    n_true = len(intensities)
    median_intensity = float(np.median(intensities)) if n_true else 1.0

    # false identifications: candidate set unrelated to the truth
    n_false = int(round(false_rate / max(1e-9, 1 - false_rate) * n_true))
    for _ in range(n_false):
        k = int(rng.integers(1, 4))
        contains = set(rng.choice(species, size=k, replace=False).tolist())
        node = tree.lca(contains)
        cand.append(frozenset(contains))
        lca_nodes.append(node)
        lca_desc.append(tree.descendant_leaves(node))
        intensities.append(median_intensity * rng.lognormal(0.0, 1.0))
        is_false.append(True)

    is_false_arr = np.array(is_false, dtype=bool)
    n_target = is_false_arr.size

    # scores: correct targets score well, false targets and decoys do not
    scores_t = np.where(
        is_false_arr,
        rng.normal(2.0, 1.0, n_target),
        rng.normal(6.0, 1.0, n_target),
    )
    n_decoy = int(rng.binomial(n_false, 0.5) + n_false * 0.5) if n_false else 0
    n_decoy = max(n_decoy, 0)
    scores_d = rng.normal(2.0, 1.0, n_decoy)

    return SimulatedRun(
        tree=tree,
        species=species,
        true_abundance=true_abundance,
        candidate_sets=cand,
        lca_nodes=lca_nodes,
        lca_descendant_sets=lca_desc,
        intensities=np.array(intensities, dtype=float),
        is_false=is_false_arr,
        scores=np.concatenate([scores_t, scores_d]),
        is_decoy=np.concatenate([np.zeros(n_target, bool), np.ones(n_decoy, bool)]),
    )
