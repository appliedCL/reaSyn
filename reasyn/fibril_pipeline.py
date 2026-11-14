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


