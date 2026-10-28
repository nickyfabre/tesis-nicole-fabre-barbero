"""Placing metaproteomic abundance on a species tree.

Modules
-------
tree      rooted tree with last-common-ancestor queries
pep       posterior error probabilities from target-decoy competition
allocate  allocation of peptide intensity to species and tree nodes (EM)
simulate  ground-truth communities for method evaluation
metrics   agreement between an estimate and a known composition
"""
from .allocate import (
    AllocationResult,
    all_rank_candidate_sets,
    allocate_species,
    node_totals_at_or_below,
    nodes_from_species,
    rigid_lca_nodes,
)
from .pep import decoy_consistency, estimate_peps
from .tree import Tree, flat_two_level_tree

__all__ = [
    "AllocationResult",
    "Tree",
    "all_rank_candidate_sets",
    "allocate_species",
    "decoy_consistency",
    "estimate_peps",
    "flat_two_level_tree",
    "node_totals_at_or_below",
    "nodes_from_species",
    "rigid_lca_nodes",
]
