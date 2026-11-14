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


def _run_analysis(root, config_path):
    root = Path(root).resolve()
    config = json.loads(Path(config_path).read_text())
    if config["causal_exposure"] is not None or config["causal_confounders"] or config["causal_target"] is not None:
        raise ValueError("causal configuration requires a validated independent-unit design; this runner supports descriptive analysis")
    for key in ("resampling_runs", "bootstrap_runs", "max_copies", "minimum_common_residues"):
        if type(config[key]) is not int or config[key] < 1:
            raise ValueError(f"{key} must be a positive integer")
    if not 0 < config["primary_cutoff"] < config["near_duplicate_score"] <= 1:
        raise ValueError("invalid label or redundancy threshold")
    if config["primary_cutoff"] not in config["label_cutoffs"]:
        raise ValueError("primary cutoff must appear in the sensitivity settings")
    manifest, curated, curation = validate_sources(root)
    rows = enrich_provenance(curated, manifest)
    code_files = sorted((root / "reasyn").glob("*.py"))
    identity = {"config": config, "source_version": manifest["dataset_version"],
                "curation_sha256": sha256_file(root / curation["path"]),
                "code": {p.name: sha256_file(p) for p in code_files},
                "environment": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "scikit-learn", "gemmi", "joblib", "threadpoolctl")},
                "python": sys.version.split()[0], "platform": platform.platform(),
                "requirements_analysis_sha256": sha256_file(root / "requirements-analysis.txt")}
    run_id = hashlib.sha256(canonical_bytes(identity)).hexdigest()[:16]
    output = root / "data/analysis" / run_id
    if (output / "completion.json").exists():
        completed = json.loads((output / "completion.json").read_text())
        for relative, checksum in completed["outputs"].items():
            if sha256_file(output / relative) != checksum:
                raise ValueError(f"cached output altered: {relative}")
        print(f"verified cached run: {output}")
        return output
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "protocol.json", identity)
    structures, features, rejected = {}, {}, []
    for row in rows:
        try:
            chains = unique_backbones(read_backbones(root / row["coordinate_path"], config["minimum_common_residues"]))
            structures[row["pdb_id"]] = chains
            features[row["pdb_id"]] = geometric_features(chains)
        except ValueError as error:
            rejected.append({"pdb_id": row["pdb_id"], "reason": str(error)})
    rows = [row for row in rows if row["pdb_id"] in structures]
    ids = [row["pdb_id"] for row in rows]
    if not set(config["reference_ids"]).issubset(ids):
        raise ValueError("a prespecified reference failed coordinate validation")
    print(f"validated {len(rows)} structures; {len(rejected)} require residue mapping", flush=True)
    write_csv(output / "coordinate_exclusions.csv", rejected, ["pdb_id", "reason"])
    write_csv(output / "backbone_inventory.csv", [
        {"pdb_id": identifier, "chain": chain["chain"], "residues": json.dumps(chain["residues"].tolist()),
         "residue_names": json.dumps(chain["names"]), "units": "angstrom", "comparison_unit": "single-chain backbone"}
        for identifier, chains in structures.items() for chain in chains])
    similarities = []
    scores = {}
    for i, left in enumerate(ids):
        for right in ids[i:]:
            result = compare_structures(structures[left], structures[right], config["minimum_common_residues"])
            similarities.append(dict(left_id=left, right_id=right, **result))
            scores[left, right] = scores[right, left] = result["score"]
        if (i + 1) % 25 == 0:
            print(f"compared {i + 1}/{len(ids)} structures", flush=True)
    write_csv(output / "similarity_matrix.csv", similarities)
    dependencies = []
    batches = curation["record"]["paper_curation"]
    if isinstance(batches, dict):
        batches = [batches]
    for batch in batches:
        for dependency in batch.get("model_dependencies", {}).values():
            if dependency.get("initial_model_pdb"):
                dependencies.extend((identifier, dependency["initial_model_pdb"]) for identifier in dependency["pdb_ids"])
    groups = dependency_groups(rows, similarities, config["near_duplicate_score"], dependencies)
    reference_groups = {groups[identifier] for identifier in config["reference_ids"]}
    group_dates = {group: min(row["release_date"] for row in rows if groups[row["pdb_id"]] == group)
                   for group in set(groups.values())}
    for row in rows:
        identifier = row["pdb_id"]
        row["dependency_group"] = groups[identifier]
        row["split"] = ("reference" if groups[identifier] in reference_groups else
                        "test" if group_dates[groups[identifier]] >= config["evaluation_from"] else "train")
        row["reference_similarity"] = max(scores[identifier, reference] for reference in config["reference_ids"])
        row["label"] = int(row["reference_similarity"] >= config["primary_cutoff"])
    write_csv(output / "frozen_splits.csv", [{key: row[key] for key in ("pdb_id", "study_id", "dependency_group", "split")} for row in rows])
    write_csv(output / "structural_metadata.csv", rows)
    write_csv(output / "geometric_features.csv", [dict(pdb_id=identifier, **values) for identifier, values in features.items()])
    train = [row for row in rows if row["split"] == "train"]
    test = [row for row in rows if row["split"] == "test"]
    if {row["dependency_group"] for row in train} & {row["dependency_group"] for row in test}:
        raise ValueError("dependency group crosses train and test")
    causal_status = {
        "eligible": False,
        "reasons": ["no agreed exposure, measured confounder set and target distribution",
                    "deposited-only data do not identify deposition probabilities",
                    "sample independence and most publication-derived conditions remain unresolved",
                    "publication annotations await human spot-check"],
        "causal_weighting": "not estimated", "causal_resampling": "not performed",
        "ordinary_study_balancing": "descriptive baseline, without a causal interpretation"}
    write_json(output / "causal_eligibility.json", causal_status)
    split_summary = {split: {"entries": len(subset), "groups": len({r["dependency_group"] for r in subset}),
                             "studies": len({r["study_id"] for r in subset}), "labels": dict(Counter(r["label"] for r in subset))}
                     for split in ("reference", "train", "test") if (subset := [r for r in rows if r["split"] == split])}
    metrics, predictions, draws_export, model_records, balance = [], [], [], [], []
    feature_names = list(next(iter(features.values())))
    for cutoff in config["label_cutoffs"]:
        y_train = np.array([int(r["reference_similarity"] >= cutoff) for r in train])
        y_test = np.array([int(r["reference_similarity"] >= cutoff) for r in test])
        if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
            metrics.append({"cutoff": cutoff, "status": "not estimable: train or test lacks both classes"})
            continue
        x_train = np.array([[features[r["pdb_id"]][name] for name in feature_names] for r in train])
        x_test = np.array([[features[r["pdb_id"]][name] for name in feature_names] for r in test])
        units = np.array([r["source_unit_id"] for r in train])
        unit_counts = Counter(units)
        group_ids = np.array([r["dependency_group"] for r in test])
        original_scaler = StandardScaler().fit(x_train)
        arms = ("original", "class_balanced", "ordinary_resampling", "study_weighted", "study_resampled", "random_weight_control")
        for arm in arms:
            repetitions = config["resampling_runs"] if arm in ("ordinary_resampling", "study_resampled", "random_weight_control") else 1
            for repetition in range(repetitions):
                seed = config["seed"] + repetition
                weights = np.ones(len(train))
                if arm == "class_balanced":
                    weights = np.array([len(train) / (2 * np.sum(y_train == value)) for value in y_train])
                elif arm == "study_weighted":
                    weights = np.array([1 / unit_counts[unit] for unit in units])
                elif arm == "random_weight_control":
                    random_weights = dict(zip(sorted(unit_counts), np.random.default_rng(seed).uniform(0.5, 1.5, len(unit_counts))))
                    weights = np.array([random_weights[unit] for unit in units])
                elif arm in ("ordinary_resampling", "study_resampled"):
                    unit_weights = {unit: count if arm == "ordinary_resampling" else 1.0 for unit, count in unit_counts.items()}
                    draws, weights = hierarchical_draw(units, unit_weights, len(train), seed, config["max_copies"])
                    if cutoff == config["primary_cutoff"]:
                        for draw in draws:
                            row = train[draw["index"]]
                            draws_export.append({"arm": arm, "seed": seed, "draw": draw["draw"],
                                                 "pdb_id": row["pdb_id"], "study_id": row["study_id"],
                                                 "sample_id": row["sample_id"], "coordinate_path": row["coordinate_path"],
                                                 "coordinate_url": row["coordinate_url"], "coordinate_sha256": row["coordinate_sha256"],
                                                 "licence": row["coordinate_licence"], "draw_count": int(weights[draw["index"]]),
                                                 "initial_probability": draw["initial_probability"],
                                                 "conditional_probability": draw["conditional_probability"]})
                weights = weights * len(train) / weights.sum()
                if cutoff == config["primary_cutoff"]:
                    for field in ("atlas_origin", "release_year", "study_id"):
                        for level in sorted({str(row[field]) for row in train}):
                            indicator = np.array([str(row[field]) == level for row in train])
                            balance.append({"arm": arm, "seed": seed, "field": field, "level": level,
                                            "original_fraction": float(indicator.mean()),
                                            "sampled_fraction": float(np.average(indicator, weights=weights)),
                                            "interpretation": "descriptive representation; no causal target specified"})
                model = LogisticRegression(C=1.0, solver="lbfgs", max_iter=2000, random_state=seed)
                model.fit(original_scaler.transform(x_train), y_train, sample_weight=weights)
                if model.n_iter_[0] >= model.max_iter:
                    raise RuntimeError("model did not converge")
                probability = model.predict_proba(original_scaler.transform(x_test))[:, 1]
                metrics.append({"cutoff": cutoff, "arm": arm, "seed": seed, "status": "exploratory",
                                "entry_weight_ess": effective_sample_size(weights), **metric_values(y_test, probability)})
                if cutoff == config["primary_cutoff"]:
                    model_records.append({"arm": arm, "seed": seed, "features": feature_names,
                                          "scaler_mean": original_scaler.mean_.tolist(), "scaler_scale": original_scaler.scale_.tolist(),
                                          "coefficients": model.coef_[0].tolist(), "intercept": float(model.intercept_[0])})
                    for row, label, value in zip(test, y_test, probability):
                        predictions.append({"arm": arm, "seed": seed, "pdb_id": row["pdb_id"],
                                            "dependency_group": row["dependency_group"], "source": row["atlas_origin"],
                                            "label": int(label), "probability": float(value)})
    write_json(output / "evaluation_metrics.json", metrics)
    write_json(output / "models.json", model_records)
    write_csv(output / "representation_diagnostics.csv", balance)
    write_csv(output / "test_predictions.csv", predictions)
    write_csv(output / "resampled_training_rows.csv", draws_export)
    intervals = {}
    for arm in sorted({r["arm"] for r in predictions}):
        arm_rows = [r for r in predictions if r["arm"] == arm]
        ordered_ids = sorted({r["pdb_id"] for r in arm_rows})
        averaged = [np.mean([r["probability"] for r in arm_rows if r["pdb_id"] == identifier]) for identifier in ordered_ids]
        labels = np.array([next(r["label"] for r in arm_rows if r["pdb_id"] == identifier) for identifier in ordered_ids])
        groups_test = np.array([next(r["dependency_group"] for r in arm_rows if r["pdb_id"] == identifier) for identifier in ordered_ids])
        intervals[arm] = bootstrap_intervals(labels, np.array(averaged), groups_test, config["seed"], config["bootstrap_runs"])
    write_json(output / "group_bootstrap.json", intervals)
    subgroup_metrics = []
    for arm in sorted({r["arm"] for r in predictions}):
        first_seed = [r for r in predictions if r["arm"] == arm and r["seed"] == config["seed"]]
        for field in ("source", "dependency_group"):
            for level in sorted({r[field] for r in first_seed}):
                subset = [r for r in first_seed if r[field] == level]
                subgroup_metrics.append({"arm": arm, "seed": config["seed"], "field": field, "level": level,
                                         "entries": len(subset), **metric_values(np.array([r["label"] for r in subset]),
                                                                                np.array([r["probability"] for r in subset]))})
    write_json(output / "subgroup_metrics.json", subgroup_metrics)
    summary = {"run_id": run_id, "included": len(curated), "analysed": len(rows), "excluded_coordinates": len(rejected),
               "splits": split_summary, "causal_eligible": False, "primary_model_fits": len(model_records),
               "total_model_fits": sum("arm" in row for row in metrics),
               "output_rows": len(draws_export), "unique_resampled_conformations": len({r["pdb_id"] for r in draws_export}),
               "evaluation_independence": "temporal and dependency-group holdout within the same atlas; external independence unverified",
               "limitations": ["single-chain comparisons omit protofilament packing and cofactors",
                               "author residue numbering is assumed homologous where no insertion flag is present",
                               "labels measure similarity to selected references, not disease biology",
                               "unrecorded shared samples, restraints and upstream fitting may remain",
                               "fixed model hyperparameters; sensitivity analyses are descriptive and do not select a winner",
                               "no demonstrated causal improvement; causal arms unavailable"]}
    write_json(output / "summary.json", summary)
    report = ["# reaSyn structural analysis", "", f"Run: `{run_id}`", "",
              f"Analysed {len(rows)} of {len(curated)} included structures. Coordinate exclusions: {len(rejected)}.", "",
              "## Evaluation", "", "Six descriptive baseline arms use identical geometric features and a fixed regularised logistic model.",
              "The temporal cohort was assigned before fitting. Reference dependency groups are excluded from model fitting and evaluation.",
              "This is an exploratory benchmark; external independence has not been established.", "",
              "| split | entries | dependency groups |", "| --- | ---: | ---: |"]
    report.extend(f"| {name} | {value['entries']} | {value['groups']} |" for name, value in split_summary.items())
    report.extend(["", "| arm | mean average precision | mean ROC AUC |", "| --- | ---: | ---: |"])
    for arm in sorted({r.get("arm") for r in metrics if "arm" in r}):
        selected = [r for r in metrics if r.get("arm") == arm and r["cutoff"] == config["primary_cutoff"]]
        if selected:
            report.append(f"| {arm} | {np.mean([r['average_precision'] for r in selected]):.3f} | {np.mean([r['roc_auc'] for r in selected]):.3f} |")
    report.extend(["", "## Causal interpretation", "", *[f"- {reason}." for reason in causal_status["reasons"]],
                   "", "No causal weights or causal dataset are released. Study balancing is a descriptive baseline.",
                   "", "## Limitations", "", *[f"- {reason}." for reason in summary["limitations"]]])
    (output / "report.md").write_text("\n".join(report) + "\n")
    write_json(output / "completion.json", {"run_id": run_id, "outputs": {
        p.relative_to(output).as_posix(): sha256_file(p) for p in sorted(output.iterdir()) if p.is_file() and p.name != "completion.json"}})
    print(json.dumps(summary, indent=2), flush=True)
    return output


def predict_saved_model(record, feature_rows):
    from scipy.special import expit
    matrix = np.array([[row[name] for name in record["features"]] for row in feature_rows], dtype=float)
    if not np.isfinite(matrix).all():
        raise ValueError("prediction features must be finite")
    scaled = (matrix - np.array(record["scaler_mean"])) / np.array(record["scaler_scale"])
    return expit(scaled @ np.array(record["coefficients"]) + record["intercept"])


