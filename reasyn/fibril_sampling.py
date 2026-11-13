"""Study-level sampling and explicit empirical back-door adjustment."""

import numpy as np


def effective_sample_size(weights):
    weights = np.asarray(weights, dtype=float)
    if weights.ndim != 1 or not np.isfinite(weights).all() or np.any(weights < 0) or weights.sum() <= 0:
        raise ValueError("weights must be finite, non-negative and have positive total")
    return float(weights.sum() ** 2 / np.dot(weights, weights))


def stratified_causal_weights(exposure, confounders, minimum_cell=2):
    """Estimate marginal-target weights for independent units with categorical C."""
    exposure = np.asarray(exposure)
    confounders = np.asarray(confounders)
    if exposure.ndim != 1 or confounders.ndim != 1 or len(exposure) != len(confounders):
        raise ValueError("one exposure and confounder stratum are required per independent unit")
    if set(exposure) != {0, 1}:
        raise ValueError("both binary exposure groups are required")
    marginal = exposure.mean()
    propensity = np.zeros(len(exposure))
    for stratum in np.unique(confounders):
        mask = confounders == stratum
        counts = np.bincount(exposure[mask].astype(int), minlength=2)
        if counts.min() < minimum_cell:
            raise ValueError(f"insufficient overlap in stratum {stratum}")
        propensity[mask] = exposure[mask].mean()
    observed = np.where(exposure == 1, propensity, 1 - propensity)
    weights = np.where(exposure == 1, marginal, 1 - marginal) / observed
    return weights, propensity


