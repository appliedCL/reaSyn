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


