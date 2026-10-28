"""Posterior error probabilities from target-decoy competition.

Estimator
---------
Let h(s) = P(decoy | search score s), estimated non-parametrically and forced
monotone. Writing N_T, N_D for the number of target and decoy PSMs and f_T, f_D
for their score densities,

    h = N_D f_D / (N_D f_D + N_T f_T)   =>   f_D / f_T = (N_T / N_D) * h / (1 - h)

With the target-decoy assumption that incorrect targets follow the decoy score
distribution, the posterior error probability is

    PEP(s) = pi0 * f_D(s) / f_T(s) = pi0 * (N_T / N_D) * h(s) / (1 - h(s))

where pi0 is the fraction of targets that are incorrect.

Under a concatenated search with target-decoy *competition* -- one best hit per
spectrum, so a spectrum contributes either a target or a decoy -- the decoy count
estimates the incorrect-target count, giving pi0 = N_D / N_T. The database-size
factors then cancel and PEP(s) = h(s) / (1 - h(s)). This is the default.

That default is what makes the estimator self-consistent with target-decoy FDR:
summing PEP over the targets in any score region returns the decoy count in that
region (`decoy_consistency`). If the search was NOT competitive -- targets and
decoys scored independently, so N_D is unrelated to the number of incorrect
targets -- pass `pi0` explicitly, or the PEPs will be inflated.
"""
from __future__ import annotations

import numpy as np


def pava_nondecreasing(y: np.ndarray, w: np.ndarray | None = None) -> np.ndarray:
    """Weighted pool-adjacent-violators: least-squares non-decreasing fit to `y`."""
    y = np.asarray(y, dtype=float)
    n = y.size
    if n == 0:
        return y.copy()
    w = np.ones(n) if w is None else np.asarray(w, dtype=float).copy()

    values = y.copy()
    weights = w.copy()
    sizes = np.ones(n, dtype=int)
    top = 0  # index of last active block

    for i in range(1, n):
        values[top + 1] = y[i]
        weights[top + 1] = w[i]
        sizes[top + 1] = 1
        top += 1
        # merge while the block order is violated
        while top > 0 and values[top - 1] > values[top]:
            tw = weights[top - 1] + weights[top]
            values[top - 1] = (weights[top - 1] * values[top - 1] + weights[top] * values[top]) / tw
            weights[top - 1] = tw
            sizes[top - 1] += sizes[top]
            top -= 1

    out = np.empty(n, dtype=float)
    pos = 0
    for b in range(top + 1):
        out[pos : pos + sizes[b]] = values[b]
        pos += sizes[b]
    return out


def estimate_peps(
    scores: np.ndarray,
    is_decoy: np.ndarray,
    *,
    pi0: float | None = None,
    higher_is_better: bool = True,
    eps: float = 1e-12,
) -> np.ndarray:
    """PEP for every PSM (targets and decoys), aligned with the inputs.

    `scores` may contain ties; PSMs with equal scores receive equal PEPs.

    `pi0` is the fraction of targets that are incorrect. Leave it None for a
    competitive target-decoy search, where it is taken as N_D / N_T; set it
    explicitly if targets and decoys were searched separately.
    """
    scores = np.asarray(scores, dtype=float)
    is_decoy = np.asarray(is_decoy).astype(bool)
    if scores.shape != is_decoy.shape:
        raise ValueError("scores and is_decoy must have the same shape")
    if scores.size == 0:
        return np.zeros(0)
    if not is_decoy.any():
        raise ValueError("no decoys present; cannot estimate PEPs by target-decoy competition")

    signed = -scores if higher_is_better else scores  # ascending = worse
    order = np.argsort(signed, kind="mergesort")
    y = is_decoy[order].astype(float)

    # collapse ties so equal scores share one fitted value
    sorted_signed = signed[order]
    boundaries = np.flatnonzero(np.diff(sorted_signed)) + 1
    group_start = np.concatenate(([0], boundaries))
    group_end = np.concatenate((boundaries, [y.size]))
    counts = (group_end - group_start).astype(float)
    means = np.add.reduceat(y, group_start) / counts

    n_t = int((~is_decoy).sum())
    n_d = int(is_decoy.sum())
    if n_t == 0:
        raise ValueError("no target PSMs present")
    scale = (n_d / n_t if pi0 is None else pi0) * (n_t / n_d)

    h_groups = pava_nondecreasing(means, counts)
    h_groups = np.clip(h_groups, 0.0, 1.0 - eps)
    pep_groups = np.clip(scale * h_groups / (1.0 - h_groups), 0.0, 1.0)

    pep_sorted = np.repeat(pep_groups, (group_end - group_start))
    out = np.empty_like(pep_sorted)
    out[order] = pep_sorted
    return out


def decoy_consistency(peps: np.ndarray, is_decoy: np.ndarray) -> tuple[float, int]:
    """(sum of target PEPs, decoy count). Equal when the estimator is consistent."""
    peps = np.asarray(peps, dtype=float)
    is_decoy = np.asarray(is_decoy).astype(bool)
    return float(peps[~is_decoy].sum()), int(is_decoy.sum())
