"""Agreement between an estimated composition and the known truth."""
from __future__ import annotations

import numpy as np


def l1_distance(est: np.ndarray, truth: np.ndarray) -> float:
    """Total variation between two compositions; 0 is perfect, 2 is disjoint."""
    return float(np.abs(np.asarray(est) - np.asarray(truth)).sum())


def mean_abs_pct_error(est: np.ndarray, truth: np.ndarray) -> float:
    """Mean absolute error in percentage points, over species present in truth."""
    est, truth = np.asarray(est), np.asarray(truth)
    present = truth > 0
    if not present.any():
        return float("nan")
    return float(np.abs(est[present] - truth[present]).mean() * 100.0)


def absent_leakage(est: np.ndarray, truth: np.ndarray) -> float:
    """Share of the estimate placed on species that are truly absent (percent)."""
    est, truth = np.asarray(est), np.asarray(truth)
    absent = truth <= 0
    return float(est[absent].sum() * 100.0) if absent.any() else 0.0


def log_pearson(est: np.ndarray, truth: np.ndarray, floor: float = 1e-9) -> float:
    est, truth = np.asarray(est), np.asarray(truth)
    keep = (truth > 0) & (est > 0)
    if keep.sum() < 3:
        return float("nan")
    a, b = np.log10(np.maximum(est[keep], floor)), np.log10(np.maximum(truth[keep], floor))
    return float(np.corrcoef(a, b)[0, 1])


def summarise(est: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    return {
        "L1": l1_distance(est, truth),
        "MAE_pct": mean_abs_pct_error(est, truth),
        "absent_leak_pct": absent_leakage(est, truth),
        "log_r": log_pearson(est, truth),
    }
