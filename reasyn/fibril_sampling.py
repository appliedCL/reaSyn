"""Study-level sampling and explicit empirical back-door adjustment."""

import numpy as np


def effective_sample_size(weights):
    weights = np.asarray(weights, dtype=float)
    if weights.ndim != 1 or not np.isfinite(weights).all() or np.any(weights < 0) or weights.sum() <= 0:
        raise ValueError("weights must be finite, non-negative and have positive total")
    return float(weights.sum() ** 2 / np.dot(weights, weights))


