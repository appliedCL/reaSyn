"""Validate frozen sources and retain field-level evidence without filling gaps."""

import csv
import hashlib
import json
import re
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    content = canonical_bytes(value)
    if path.exists() and path.read_bytes() == content:
        return
    temporary = path.with_suffix(path.suffix + ".part")
    temporary.write_bytes(content)
    temporary.replace(path)


def write_csv(path, rows, fields=None):
    rows = list(rows)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else [])
    temporary = path.with_suffix(path.suffix + ".part")
    with temporary.open("w", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def read_csv(path):
    with Path(path).open(newline="") as source:
        return list(csv.DictReader(source))


def project_path(root, relative):
    root = Path(root).resolve()
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f"source path escapes the project: {relative}")
    return path


def validate_sources(root):
    root = Path(root)
    manifest = json.loads((root / "data/fibril_manifest.json").read_text())
    if manifest.get("failures"):
        raise ValueError("source manifest records unresolved download failures")
    for relative, record in manifest["files"].items():
        path = project_path(root, relative)
        if not path.is_file() or sha256_file(path) != record["sha256"]:
            raise ValueError(f"missing or altered source: {relative}")
    computed_version = hashlib.sha256("\n".join(
        f"{name}:{record['sha256']}" for name, record in sorted(manifest["files"].items())
    ).encode()).hexdigest()
    if computed_version != manifest["dataset_version"]:
        raise ValueError("source dataset version does not match its checksums")
    review = read_csv(root / "data/fibril_inclusion_review.csv")
    if len({row["pdb_id"] for row in review}) != len(review):
        raise ValueError("duplicate inclusion-review identifiers")
    included = sorted(row["pdb_id"] for row in review if row["decision"] == "include")
    if not included:
        raise ValueError("no structures approved by the inclusion review")
    notebook = json.loads((root / "notebooks/curate_fibril_metadata.ipynb").read_text())
    code = "\n\n".join("".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code")
    code_hash = hashlib.sha256(code.encode()).hexdigest()
    matches = []
    for path in (root / "data/curated").glob("fibril_metadata_*.json"):
        record = json.loads(path.read_text())
        if record.get("notebook_code_sha256") == code_hash:
            matches.append((path, record))
    if len(matches) != 1:
        raise ValueError("run curate_fibril_metadata.ipynb: exactly one export must match its current code")
    path, record = matches[0]
    checks = {
        "source_manifest_sha256": root / "data/fibril_manifest.json",
        "inclusion_review_sha256": root / "data/fibril_inclusion_review.csv",
        "requirements_sha256": root / "requirements.txt",
        "table_sha256": path.with_suffix(".csv"),
    }
    for key, source in checks.items():
        if sha256_file(source) != record[key]:
            raise ValueError(f"stale curation export: {key}; rerun the curation notebook")
    rows = read_csv(path.with_suffix(".csv"))
    if sorted(row["pdb_id"] for row in rows) != included or len(rows) != record["rows"]:
        raise ValueError("curated rows do not match the inclusion review")
    batches = record.get("paper_curation", [])
    if isinstance(batches, dict):
        batches = [batches]
    for batch in batches:
        for source in batch.get("sources", {}).values():
            if sha256_file(project_path(root, source["path"])) != source["sha256"]:
                raise ValueError(f"altered publication: {source['path']}")
    return manifest, rows, {"path": path.relative_to(root).as_posix(), "record": record}


