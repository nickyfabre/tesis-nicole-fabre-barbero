"""Allocating peptide intensity to species, and from species to tree nodes.

Model
-----
Within one run, each species s carries a share pi_s of the community's
detectable protein, sum_s pi_s = 1. Observation i has intensity I_i, a candidate
set C_i of species whose proteomes contain the peptide, and a posterior error
probability p_i that the identification is wrong. Treating intensity as a weight
and letting a correct observation arise from exactly one candidate species,

    L(pi) = sum_i I_i (1 - p_i) log( sum_{s in C_i} pi_s )

The EM update for this weighted multinomial mixture is

    pi_s  <-  pi_s * sum_{C containing s} W_C / (sum_{s' in C} pi_s')      (then renormalise)
    W_C   =   sum_{i : C_i = C} I_i (1 - p_i)

which is the estimator used for multi-mapping reads in RSEM and Salmon, with
reads replaced by intensity and downweighted by identification confidence.

The strategies compared in the paper are nested special cases:

  unique        only observations with |C_i| = 1 contribute
  proportional  one EM iteration started from the unique-peptide estimate
  em            the same iteration run to convergence
  em+PEP        as em, with p_i from target-decoy competition (pass `peps`)

Because the likelihood depends on the observations only through W_C, the
iteration runs over distinct candidate sets rather than over observations.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

from .tree import Tree


@dataclass
class AllocationResult:
    species: list[str]
    pi: np.ndarray
    n_iter: int
    converged: bool
    loglik: float
    loglik_trace: list[float] = field(default_factory=list)

    def as_dict(self) -> dict[str, float]:
        return {s: float(v) for s, v in zip(self.species, self.pi)}


def _group_candidate_sets(
    candidate_sets: Sequence[frozenset[str]],
    intensities: np.ndarray,
    species_index: Mapping[str, int],
    peps: np.ndarray | None,
) -> tuple[list[np.ndarray], np.ndarray]:
    """Collapse observations to (distinct candidate set, total weight)."""
    intensities = np.asarray(intensities, dtype=float)
    if peps is None:
        weights = intensities
    else:
        peps = np.asarray(peps, dtype=float)
        if peps.shape != intensities.shape:
            raise ValueError("peps and intensities must have the same shape")
        weights = intensities * (1.0 - peps)

    totals: dict[frozenset[str], float] = {}
    for cset, w in zip(candidate_sets, weights):
        if w <= 0.0 or not cset:
            continue
        known = cset & species_index.keys()
        if not known:
            continue
        key = frozenset(known)
        totals[key] = totals.get(key, 0.0) + float(w)

    sets = [np.array(sorted(species_index[s] for s in cset), dtype=np.intp) for cset in totals]
    return sets, np.array(list(totals.values()), dtype=float)


def _loglik(sets: Sequence[np.ndarray], weights: np.ndarray, pi: np.ndarray) -> float:
    total = 0.0
    for idx, w in zip(sets, weights):
        total += w * np.log(max(pi[idx].sum(), 1e-300))
    return float(total)


def allocate_species(
    candidate_sets: Sequence[frozenset[str]],
    intensities: np.ndarray,
    species: Sequence[str],
    *,
    method: str = "em",
    peps: np.ndarray | None = None,
    init: str = "uniform",
    max_iter: int = 1000,
    tol: float = 1e-10,
    floor: float = 1e-12,
) -> AllocationResult:
    """Estimate per-species shares within one run.

    `method` is one of 'unique', 'proportional', 'em'. Passing `peps` weights
    each observation by its probability of being a correct identification.
    `init` is 'uniform' or 'unique'; with 'unique' a species having no unique
    peptide starts at the floor and can still gain signal, which a single
    proportional step cannot give it.
    """
    if method not in {"unique", "proportional", "em"}:
        raise ValueError(f"unknown method {method!r}")
    species = list(species)
    index = {s: i for i, s in enumerate(species)}
    n = len(species)
    sets, weights = _group_candidate_sets(candidate_sets, intensities, index, peps)

    if not sets:
        return AllocationResult(species, np.zeros(n), 0, True, 0.0, [])

    # unique-peptide estimate: the anchor, and the 'unique' arm's answer
    anchors = np.zeros(n)
    for idx, w in zip(sets, weights):
        if idx.size == 1:
            anchors[idx[0]] += w

    if method == "unique":
        total = anchors.sum()
        pi = anchors / total if total > 0 else np.full(n, 1.0 / n)
        return AllocationResult(species, pi, 0, True, _loglik(sets, weights, np.maximum(pi, floor)), [])

    if init == "unique" or method == "proportional":
        pi = anchors.copy()
        if pi.sum() <= 0:
            pi = np.full(n, 1.0 / n)
        pi = pi / pi.sum()
        pi = np.maximum(pi, floor)
        pi /= pi.sum()
    elif init == "uniform":
        pi = np.full(n, 1.0 / n)
    else:
        raise ValueError(f"unknown init {init!r}")

    n_steps = 1 if method == "proportional" else max_iter
    trace = [_loglik(sets, weights, pi)]
    converged = False
    it = 0
    for it in range(1, n_steps + 1):
        numer = np.zeros(n)
        for idx, w in zip(sets, weights):
            denom = pi[idx].sum()
            if denom <= 0:
                continue
            numer[idx] += w * pi[idx] / denom
        total = numer.sum()
        if total <= 0:
            break
        new_pi = numer / total
        shift = float(np.abs(new_pi - pi).sum())
        pi = new_pi
        trace.append(_loglik(sets, weights, np.maximum(pi, 1e-300)))
        if method == "em" and shift < tol:
            converged = True
            break
    if method == "proportional":
        converged = True

    return AllocationResult(species, pi, it, converged, trace[-1], trace)


def nodes_from_species(
    tree: Tree,
    species_share: Mapping[str, float],
    nodes: Iterable[str] | None = None,
) -> dict[str, float]:
    """Roll leaf-level shares up to every node (a node = sum over its leaves)."""
    targets = list(tree.nodes()) if nodes is None else list(nodes)
    return {
        node: float(sum(species_share.get(leaf, 0.0) for leaf in tree.descendant_leaves(node)))
        for node in targets
    }


def rigid_lca_nodes(
    tree: Tree,
    assigned_node: Sequence[str],
    intensities: np.ndarray,
    *,
    peps: np.ndarray | None = None,
) -> dict[str, float]:
    """Retain-at-LCA: a node's quantity is the signal assigned at or below it.

    No redistribution, so this is the assumption-free view of the same evidence.
    """
    intensities = np.asarray(intensities, dtype=float)
    weights = intensities if peps is None else intensities * (1.0 - np.asarray(peps, dtype=float))

    direct: dict[str, float] = {}
    for node, w in zip(assigned_node, weights):
        if w > 0:
            direct[node] = direct.get(node, 0.0) + float(w)

    out = {node: 0.0 for node in tree.nodes()}
    for node, w in direct.items():
        for ancestor in tree.ancestors(node):
            out[ancestor] += w
    total = out[tree.root]
    if total > 0:
        out = {k: v / total for k, v in out.items()}
    return out


def all_rank_candidate_sets(
    tree: Tree, candidate_sets: Sequence[frozenset[str]]
) -> list[frozenset[str]]:
    """Candidate *nodes* for an all-ranks model, not just leaves.

    Reimplements the modelling choice made by MiCId, which estimates abundance
    for taxa at every rank rather than placing all mass on leaves: a peptide
    shared across a genus is evidence for the genus taxon in its own right as
    well as for each species in it. For a peptide contained in species set C the
    candidate nodes are those lying on a path from LCA(C) down to some member of
    C. This is our reimplementation of that idea, not the MiCId software.
    """
    out: list[frozenset[str]] = []
    cache: dict[frozenset[str], frozenset[str]] = {}
    for cset in candidate_sets:
        if cset in cache:
            out.append(cache[cset])
            continue
        root = tree.lca(cset)
        nodes: set[str] = set()
        for leaf in cset:
            for anc in tree.ancestors(leaf):
                nodes.add(anc)
                if anc == root:
                    break
        result = frozenset(nodes)
        cache[cset] = result
        out.append(result)
    return out


def node_totals_at_or_below(tree: Tree, node_share: Mapping[str, float]) -> dict[str, float]:
    """Aggregate an all-ranks estimate: a node gets its own mass plus everything below."""
    out = {n: 0.0 for n in tree.nodes()}
    for node, value in node_share.items():
        if value:
            for ancestor in tree.ancestors(node):
                out[ancestor] += float(value)
    return out
