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


def hierarchical_draw(units, unit_weights, size, seed, max_copies=5, within_weights=None):
    """Return sequential capped draws and their conditional selection probabilities."""
    units = np.asarray(units)
    if size < 1 or max_copies < 1 or size > len(units) * max_copies:
        raise ValueError("requested draw cannot satisfy the copy cap")
    unique = np.unique(units)
    weights = np.array([unit_weights[unit] for unit in unique], dtype=float)
    effective_sample_size(weights)
    population = np.ones(len(units)) if within_weights is None else np.asarray(within_weights, dtype=float)
    if population.shape != units.shape or not np.isfinite(population).all() or np.any(population < 0):
        raise ValueError("invalid within-unit population weights")
    indices = [np.flatnonzero(units == unit) for unit in unique]
    if any(population[idx].sum() <= 0 for idx in indices):
        raise ValueError("each unit needs positive within-unit mass")
    capacity = sum(max_copies * np.count_nonzero(population[idx]) for idx, weight in zip(indices, weights) if weight > 0)
    if size > capacity:
        raise ValueError("positive-weight conformations cannot satisfy the requested cap")
    rng = np.random.default_rng(seed)
    counts = np.zeros(len(units), dtype=int)
    draws = []
    initial = np.zeros(len(units))
    for idx, weight in zip(indices, weights):
        initial[idx] = weight / weights.sum() * population[idx] / population[idx].sum()
    for draw in range(size):
        available = [(counts[idx] < max_copies) & (population[idx] > 0) for idx in indices]
        active = weights * np.array([mask.any() for mask in available])
        active /= active.sum()
        unit_index = rng.choice(len(unique), p=active)
        idx = indices[unit_index]
        conditional = population[idx] * available[unit_index]
        conditional /= conditional.sum()
        position = rng.choice(len(idx), p=conditional)
        selected = idx[position]
        counts[selected] += 1
        draws.append({"draw": draw, "index": int(selected), "initial_probability": float(initial[selected]),
                      "conditional_probability": float(active[unit_index] * conditional[position])})
    return draws, counts
