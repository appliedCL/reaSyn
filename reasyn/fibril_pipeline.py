"""Offline, versioned structural analysis with explicit scientific eligibility gates."""

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import os
import platform
from pathlib import Path
import sys

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.preprocessing import StandardScaler

from .fibril_geometry import read_backbones, unique_backbones, compare_structures, geometric_features
from .fibril_sampling import hierarchical_draw, effective_sample_size
from .fibril_sources import canonical_bytes, enrich_provenance, read_csv, sha256_file, validate_sources, write_csv, write_json


def dependency_groups(rows, similarities, threshold, dependencies=()):
    parent = {row["pdb_id"]: row["pdb_id"] for row in rows}

    def find(item):
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def join(left, right):
        if left in parent and right in parent:
            a, b = sorted((find(left), find(right)))
            parent[b] = a

    studies = {}
    for row in rows:
        previous = studies.setdefault(row["study_id"], row["pdb_id"])
        join(previous, row["pdb_id"])
    for pair in similarities:
        if pair["score"] >= threshold:
            join(pair["left_id"], pair["right_id"])
    for left, right in dependencies:
        join(left, right)
    return {identifier: find(identifier) for identifier in parent}


def metric_values(labels, probabilities, weights=None):
    values = {"brier": float(brier_score_loss(labels, probabilities, sample_weight=weights))}
    if len(np.unique(labels)) == 2:
        values.update(average_precision=float(average_precision_score(labels, probabilities, sample_weight=weights)),
                      roc_auc=float(roc_auc_score(labels, probabilities, sample_weight=weights)))
    else:
        values.update(average_precision=None, roc_auc=None)
    return values


def bootstrap_intervals(labels, probabilities, groups, seed, repetitions):
    rng = np.random.default_rng(seed)
    unique = np.unique(groups)
    values = []
    for _ in range(repetitions):
        sampled = rng.choice(unique, len(unique), replace=True)
        indices = np.concatenate([np.flatnonzero(groups == group) for group in sampled])
        if len(np.unique(labels[indices])) == 2:
            values.append(metric_values(labels[indices], probabilities[indices]))
    return {"valid_replicates": len(values), "independent_groups": len(unique),
            "intervals": {metric: np.quantile([value[metric] for value in values], [0.025, 0.975]).tolist()
                          if len(unique) >= 3 and len(values) >= 30 else None
                          for metric in ("average_precision", "roc_auc", "brier")}}


def run_analysis(root, config_path):
    root = Path(root).resolve()
    lock = root / "data/.fibril-analysis.lock"
    if not lock.parent.is_dir():
        raise ValueError("download and curate the source dataset before running analysis")
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error:
        raise RuntimeError("analysis already running; remove a stale data/.fibril-analysis.lock only after checking its process") from error
    try:
        with os.fdopen(descriptor, "w") as destination:
            destination.write(str(os.getpid()))
        return _run_analysis(root, config_path)
    finally:
        lock.unlink(missing_ok=True)


