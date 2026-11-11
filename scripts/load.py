import argparse
import json
import struct
import tempfile
import zipfile
from pathlib import Path

import MDAnalysis as mda
import numpy as np

from inventory import ROOT, verify_archive


def copy_frames(source, destination, count):
    """Copy complete compressed XTC frames without changing their coordinates."""
    for _ in range(count):
        header = source.read(92)
        if len(header) != 92 or struct.unpack_from(">i", header)[0] != 1995:
            raise ValueError("Missing or invalid XTC frame header.")
        length = struct.unpack_from(">i", header, 88)[0]
        if length <= 0:
            raise ValueError("Invalid XTC payload length.")
        destination.write(header)
        remaining = (length + 3) & ~3
        while remaining:
            chunk = source.read(min(remaining, 1024 * 1024))
            if not chunk:
                raise ValueError("Truncated XTC frame.")
            destination.write(chunk)
            remaining -= len(chunk)


def load_frames(inventory_path, count=3):
    """Yield coordinate copies for the first count frames of each trajectory."""
    if count < 1:
        raise ValueError("The number of frames must be positive.")
    inventory = json.loads(Path(inventory_path).read_text())
    units = inventory["trajectories"]
    if not units or any(count > unit["trajectory"]["frame_count"] for unit in units):
        raise ValueError("Requested frames exceed the available trajectories.")
    source = inventory["archive"]
    archive_path = ROOT / source["path"]
    verify_archive(archive_path, {
        "size_bytes": source["size_bytes"], "checksum": source["published_checksum"],
    })
    with zipfile.ZipFile(archive_path) as archive:
        for unit in units:
            with tempfile.TemporaryDirectory(prefix="reasyn-load-") as directory:
                topology = Path(directory) / "reference.gro"
                trajectory = Path(directory) / "frames.xtc"
                topology.write_bytes(archive.read(unit["files"]["gro"]))
                with archive.open(unit["files"]["xtc"]) as stream, trajectory.open("wb") as target:
                    copy_frames(stream, target, count)
                universe = mda.Universe(str(topology), convert_units=False, to_guess=())
                atoms = unit["observed"]["atom_count"]
                if universe.atoms.n_atoms != atoms:
                    raise ValueError(f"{unit['trajectory_id']}: topology atom count differs.")
                reference = universe.atoms.positions.copy()
                universe.load_new(str(trajectory), convert_units=False)
                with universe.trajectory as reader:
                    if len(reader) != count:
                        raise ValueError("Decoded frame count differs from the requested count.")
                    for timestep in reader:
                        positions = timestep.positions.copy()
                        box = timestep.triclinic_dimensions
                        if positions.shape != (atoms, 3) or not np.isfinite(positions).all():
                            raise ValueError("Invalid coordinate shape or non-finite coordinates.")
                        if box is None or not np.isfinite(box).all():
                            raise ValueError("Missing or invalid periodic box.")
                        if timestep.frame == 0 and not np.allclose(
                            positions, reference, rtol=0, atol=0.001
                        ):
                            raise ValueError("First trajectory frame and GRO reference disagree.")
                        yield {
                            "structure_id": f"{unit['trajectory_id']}::frame:{timestep.frame}",
                            "source_accession": inventory["source_study"]["accession"],
                            "study_publication": inventory["source_study"]["publication"],
                            "trajectory_id": unit["trajectory_id"],
                            "condition_label": unit["condition_label"],
                            "original_sample_id": unit["unresolved_metadata"]["original_sample_id"],
                            "original_run_id": unit["unresolved_metadata"]["original_run_id"],
                            "original_replica_id": unit["unresolved_metadata"]["original_replica_id"],
                            "frame_index": timestep.frame,
                            "time_ps": float(timestep.time),
                            "archive_path": source["path"],
                            "topology_member": unit["files"]["gro"],
                            "trajectory_member": unit["files"]["xtc"],
                            "coordinates_nm": positions,
                            "box_vectors_nm": box.copy(),
                        }


def main():
    parser = argparse.ArgumentParser(description="Load a few development conformations in nm.")
    parser.add_argument("--inventory", type=Path, default=ROOT / "data/development_inventory.json")
    parser.add_argument("--frames", type=int, default=3, help="First frames to load per trajectory.")
    args = parser.parse_args()
    try:
        for frame in load_frames(args.inventory, args.frames):
            print(f"{frame['structure_id']} | {frame['time_ps']:g} ps | "
                  f"{len(frame['coordinates_nm'])} atoms | nm")
    except (OSError, ValueError, KeyError, struct.error, zipfile.BadZipFile) as error:
        parser.exit(1, f"Loading failed: {error}\n")


if __name__ == "__main__":
    main()
