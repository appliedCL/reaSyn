import argparse
import hashlib
import json
import math
import struct
import zipfile
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).resolve().parents[1]
REFERENCE_SEQUENCE = (
    "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK"
    "EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP"
    "DNEAYEMPSEEGYQDYEPEA"
)
RESIDUE_CODES = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV",
))
RESIDUE_CODES.update(HIE="H", SEP="S")


def verify_archive(path, source):
    if path.stat().st_size != source["size_bytes"]:
        raise ValueError(f"Archive size mismatch: {path}")
    md5, sha256 = hashlib.md5(), hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            md5.update(chunk)
            sha256.update(chunk)
    if "md5:" + md5.hexdigest() != source["checksum"]:
        raise ValueError(f"Archive checksum mismatch: {path}")
    return {"size_bytes": path.stat().st_size, "md5": md5.hexdigest(),
            "sha256": sha256.hexdigest()}


def reference_summary(pdb, gro, label):
    atoms = [line for line in pdb.decode("ascii").splitlines()
             if line.startswith(("ATOM  ", "HETATM"))]
    residues = [line for line in atoms if line[12:16].strip() == "CA"]
    if [int(line[22:26]) for line in residues] != list(range(1, 141)):
        raise ValueError(f"{label}: expected protein residues numbered 1–140.")
    if len({line[21:22] for line in residues}) != 1:
        raise ValueError(f"{label}: expected one protein chain.")
    sequence = "".join(RESIDUE_CODES.get(line[17:20].strip(), "?") for line in residues)
    if sequence != REFERENCE_SEQUENCE:
        raise ValueError(f"{label}: normalised sequence differs from alpha-synuclein P37840.")
    gro_lines = gro.decode("ascii").splitlines()
    if len(gro_lines) < 3 or int(gro_lines[1]) != len(atoms):
        raise ValueError(f"{label}: PDB/GRO atom counts differ.")
    return {
        "atom_count": len(atoms),
        "protein_residue_count": len(residues),
        "sequence": sequence,
        "sequence_normalisation": {"HIE": "H", "SEP": "S"},
        "s129_phosphorylated": residues[128][17:20].strip() == "SEP",
        "acetyl_cap_present": any(line[17:20].strip() == "ACE" for line in atoms),
    }


def trajectory_summary(stream, expected_atoms, label):
    """Read compressed XTC frame headers without decoding coordinates."""
    frames = 0
    first_time = last_time = minimum_interval = maximum_interval = None
    while True:
        header = stream.read(92)
        if not header:
            break
        if len(header) != 92:
            raise ValueError(f"{label}: truncated XTC header.")
        magic, atoms, _, time_ps = struct.unpack_from(">iiif", header)
        coordinate_count = struct.unpack_from(">i", header, 52)[0]
        if magic != 1995:
            raise ValueError(f"{label}: invalid XTC magic.")
        if atoms != expected_atoms or coordinate_count != atoms:
            raise ValueError(f"{label}: XTC/reference atom counts differ.")
        if not math.isfinite(time_ps):
            raise ValueError(f"{label}: non-finite XTC time.")
        length = struct.unpack_from(">i", header, 88)[0]
        if length <= 0:
            raise ValueError(f"{label}: invalid XTC payload length.")
        remaining = (length + 3) & ~3
        while remaining:
            chunk = stream.read(min(remaining, 1024 * 1024))
            if not chunk:
                raise ValueError(f"{label}: truncated XTC payload.")
            remaining -= len(chunk)
        if last_time is None:
            first_time = time_ps
        else:
            interval = time_ps - last_time
            if interval <= 0:
                raise ValueError(f"{label}: non-increasing frame time.")
            minimum_interval = interval if minimum_interval is None else min(minimum_interval, interval)
            maximum_interval = interval if maximum_interval is None else max(maximum_interval, interval)
        last_time = time_ps
        frames += 1
    if not frames:
        raise ValueError(f"{label}: empty XTC trajectory.")
    return {
        "frame_count": frames, "first_time_ps": first_time, "last_time_ps": last_time,
        "minimum_frame_interval_ps": minimum_interval,
        "maximum_frame_interval_ps": maximum_interval,
    }


def create_inventory(path, dataset, record, source):
    checksums = verify_archive(path, source)
    trajectories = []
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Archive contains duplicate member names.")
        for member in sorted(name for name in names if name.endswith(".xtc")):
            files = {ext: member[:-4] + "." + ext for ext in ("pdb", "gro", "tpr", "xtc")}
            missing = set(files.values()) - set(names)
            if missing:
                raise ValueError(f"Missing matching reference files: {sorted(missing)}")
            observed = reference_summary(
                archive.read(files["pdb"]), archive.read(files["gro"]), member
            )
            with archive.open(member) as stream:
                timing = trajectory_summary(stream, observed["atom_count"], member)
            folder = PurePosixPath(member).parent
            states = {"AS": False, "pAS": True}
            if folder.name in states and observed["s129_phosphorylated"] != states[folder.name]:
                raise ValueError(f"{member}: condition label and phosphorylation state disagree.")
            trajectories.append({
                "trajectory_id": record["accession"] + "::" + member,
                "condition_label": folder.name,
                "force_field_folder_label": folder.parts[0].removesuffix(" simulations"),
                "files": files, "observed": observed, "trajectory": timing,
                "unresolved_metadata": {
                    "original_sample_id": None, "original_run_id": None,
                    "original_replica_id": None, "pH": None, "run_seed": None,
                    "batch_id": None, "original_population_weights": None,
                },
            })
    if not trajectories:
        raise ValueError("Archive contains no XTC trajectories.")
    return {
        "schema_version": 1,
        "source_study": {
            "dataset_id": dataset["id"], "publication": dataset["publication"],
            "accession": record["accession"], "record_url": record["record_url"],
            "licence": record["licence"],
        },
        "archive": {"path": source["path"], "url": source["url"],
                    "published_checksum": source["checksum"], **checksums},
        "reported_study_conditions": {
            "temperature_k": 300, "nominal_nacl_mM": 150, "sampling_method": "REST2",
            "evidence": [record["record_url"], dataset["publication"]],
            "verified_from_simulation_parameters": False,
        },
        "trajectory_count": len(trajectories),
        "total_stored_frames": sum(item["trajectory"]["frame_count"] for item in trajectories),
        "trajectories": trajectories,
        "limitations": [
            "Stored frames are not established as unique conformations or independent samples.",
            "Coordinate geometry, PDB/GRO atom order and TPR parameters remain unchecked.",
            "A causal model and suitable confounders have not been established.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description="Create an inventory of the development archive.")
    parser.add_argument("--manifest", type=Path, default=ROOT / "data/sources.json")
    parser.add_argument("--output", type=Path, default=ROOT / "data/development_inventory.json")
    args = parser.parse_args()
    try:
        manifest_bytes = args.manifest.read_bytes()
        manifest = json.loads(manifest_bytes)
        datasets = [item for item in manifest["datasets"] if item["id"] == "de_bruyn_2024"]
        if len(datasets) != 1:
            raise ValueError("Expected one de_bruyn_2024 dataset in the manifest.")
        dataset = datasets[0]
        sources = [(record, item) for record in dataset["records"] for item in record["files"]
                   if item["selected"] and item["name"] == "protein+ion_trajectories.zip"]
        if len(sources) != 1:
            raise ValueError("Expected one selected protein-and-ion archive.")
        record, source = sources[0]
        path = ROOT / source["path"]
        if args.output.resolve() in (path.resolve(), args.manifest.resolve()):
            raise ValueError("The output must not overwrite the source archive or manifest.")
        report = create_inventory(path, dataset, record, source)
        report["source_manifest_sha256"] = hashlib.sha256(manifest_bytes).hexdigest()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(args.output.name + ".tmp")
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        temporary.replace(args.output)
        print(f"{report['trajectory_count']} trajectories; {report['total_stored_frames']:,} frames.")
        print(f"Report: {args.output}")
    except (OSError, ValueError, KeyError, TypeError, struct.error, zipfile.BadZipFile) as error:
        parser.exit(1, f"Inventory failed: {error}\n")


if __name__ == "__main__":
    main()
