import numpy as np
import pytest
import json
from pathlib import Path

from reasyn.fibril_geometry import aligned_distances, compare_backbones
from reasyn.fibril_pipeline import dependency_groups
from reasyn.fibril_sampling import hierarchical_draw, stratified_causal_weights
from reasyn.fibril_sources import classify_origin, project_path
from reasyn.fibril_sources import read_csv
from reasyn.fibril_pipeline import run_analysis, predict_saved_model


def test_rigid_transform_invariance_and_reflection():
    rng = np.random.default_rng(10)
    xyz = rng.normal(size=(40, 3))
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    rotation[:, -1] *= np.linalg.det(rotation)
    moved = xyz @ rotation + 12
    assert np.max(aligned_distances(xyz, moved)) < 1e-12
    assert np.mean(aligned_distances(xyz, -xyz)) > 0.1
    left = {"xyz": xyz, "residues": np.arange(40)}
    right = {"xyz": moved, "residues": np.arange(40)}
    assert compare_backbones(left, right)["score"] == pytest.approx(1)


