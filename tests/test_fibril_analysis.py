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


def test_missing_residues_cannot_score_as_complete():
    xyz = np.random.default_rng(4).normal(size=(60, 3))
    left = {"xyz": xyz, "residues": np.arange(60)}
    right = {"xyz": xyz[:30], "residues": np.arange(30)}
    assert compare_backbones(left, right)["score"] == pytest.approx(0.5)


def test_causal_adjustment_recovers_planted_confounding():
    c = np.repeat([0, 1], 1000)
    a = np.r_[np.ones(100), np.zeros(900), np.ones(900), np.zeros(100)].astype(int)
    outcome = 2 * a + 4 * c
    weights, propensity = stratified_causal_weights(a, c)
    naive = outcome[a == 1].mean() - outcome[a == 0].mean()
    corrected = np.average(outcome[a == 1], weights=weights[a == 1]) - np.average(outcome[a == 0], weights=weights[a == 0])
    assert naive > 5
    assert corrected == pytest.approx(2)
    assert set(np.round(propensity, 1)) == {0.1, 0.9}


def test_no_overlap_is_rejected():
    with pytest.raises(ValueError, match="overlap"):
        stratified_causal_weights([0, 0, 1, 1], [0, 0, 1, 1])


def test_constant_weights_when_exposure_is_balanced():
    weights, _ = stratified_causal_weights([0, 0, 1, 1] * 2, [0] * 4 + [1] * 4)
    assert np.all(weights == 1)


@pytest.mark.parametrize("seed", range(20))
def test_capped_sampling_determinism_and_probabilities(seed):
    units = np.array(["a", "a", "b"])
    draws, counts = hierarchical_draw(units, {"a": 1, "b": 1}, 6, seed, max_copies=2)
    again, again_counts = hierarchical_draw(units, {"a": 1, "b": 1}, 6, seed, max_copies=2)
    assert draws == again
    assert np.array_equal(counts, again_counts)
    assert counts.sum() == 6 and counts.max() == 2
    assert all(0 < d["conditional_probability"] <= 1 for d in draws)
    assert all(d["initial_probability"] == (0.5 if d["index"] == 2 else 0.25) for d in draws)


def test_impossible_cap_and_zero_mass_rejected():
    with pytest.raises(ValueError):
        hierarchical_draw(["a"], {"a": 1}, 3, 1, 2)
    with pytest.raises(ValueError):
        hierarchical_draw(["a", "b"], {"a": 1, "b": 0}, 4, 1, 2)


def test_transitive_dependencies_cannot_cross_splits():
    rows = [{"pdb_id": "a", "study_id": "one"}, {"pdb_id": "b", "study_id": "one"},
            {"pdb_id": "c", "study_id": "two"}, {"pdb_id": "d", "study_id": "three"}]
    groups = dependency_groups(rows, [{"left_id": "b", "right_id": "c", "score": 0.99}], 0.95, [("c", "d")])
    assert len(set(groups.values())) == 1


def test_clinical_seed_does_not_make_recombinant_material_patient_derived():
    assert classify_origin("human recombinant seeded from patient brain")[0] == "in_vitro"
    assert classify_origin("human extracted from M83 mouse brain")[0] == "model_organism"
    assert classify_origin("human extracted from patient brain")[0] == "human_tissue"
    assert classify_origin("human")[0] == "unknown"


def test_manifest_cannot_escape_project(tmp_path):
    with pytest.raises(ValueError, match="escapes"):
        project_path(tmp_path, "../outside")


