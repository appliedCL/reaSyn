"""Residue-matched backbone comparisons; scores are not optimised TM-scores."""

import re

import gemmi
import numpy as np
from scipy.spatial.distance import pdist, squareform


def read_backbones(path, minimum_residues=30):
    block = gemmi.cif.read_file(str(path)).sole_block()
    entities = block.get_mmcif_category("_entity.")
    selected = {identifier for identifier, description in zip(entities["id"], entities["pdbx_description"])
                if re.search(r"(alpha|α)[ -]?synuclein", description, re.I)}
    atoms = block.get_mmcif_category("_atom_site.")
    chains = {}
    for i, atom in enumerate(atoms["label_atom_id"]):
        if atom != "CA" or atoms["label_entity_id"][i] not in selected:
            continue
        if atoms["pdbx_PDB_model_num"][i] != "1":
            continue
        if atoms["pdbx_PDB_ins_code"][i] not in (None, False, ".", "?"):
            raise ValueError("insertion codes require an explicit sequence mapping")
        residue = int(atoms["auth_seq_id"][i])
        if not 1 <= residue <= 140:
            raise ValueError("author numbering outside human 1–140 requires mapping")
        key = atoms["label_asym_id"][i]
        coordinates = np.array([float(atoms[field][i]) for field in ("Cartn_x", "Cartn_y", "Cartn_z")])
        occupancy = float(atoms["occupancy"][i])
        if not np.isfinite(coordinates).all() or not np.isfinite(occupancy):
            raise ValueError("non-finite coordinates or occupancy")
        if occupancy <= 0:
            continue
        alt = str(atoms["label_alt_id"][i] or "")
        candidate = (occupancy, alt, coordinates, atoms["label_comp_id"][i])
        previous = chains.setdefault(key, {}).get(residue)
        if previous is None or (-occupancy, alt) < (-previous[0], previous[1]):
            chains[key][residue] = candidate
    result = []
    for chain, residues in sorted(chains.items()):
        if len(residues) < minimum_residues:
            continue
        indices = np.array(sorted(residues))
        xyz = np.array([residues[i][2] for i in indices])
        names = [residues[i][3] for i in indices]
        result.append({"chain": chain, "residues": indices, "xyz": xyz, "names": names})
    if not result:
        raise ValueError("no sufficiently long alpha-synuclein backbone")
    return result


def aligned_distances(left, right):
    if left.shape != right.shape or left.ndim != 2 or left.shape[1] != 3 or len(left) < 3:
        raise ValueError("alignment needs matching arrays of at least three xyz positions")
    left = left - left.mean(axis=0)
    right = right - right.mean(axis=0)
    u, _, vt = np.linalg.svd(left.T @ right)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    return np.linalg.norm(left @ (u @ correction @ vt) - right, axis=1)


def compare_backbones(left, right, minimum_common=30):
    common, li, ri = np.intersect1d(left["residues"], right["residues"], return_indices=True)
    if len(common) < minimum_common:
        return {"score": 0.0, "rmsd": None, "common": len(common), "coverage": 0.0}
    distances = aligned_distances(left["xyz"][li], right["xyz"][ri])
    length = max(len(left["xyz"]), len(right["xyz"]))
    scale = max(0.5, 1.24 * np.cbrt(max(length - 15, 1)) - 1.8)
    return {"score": float(np.sum(1 / (1 + (distances / scale) ** 2)) / length),
            "rmsd": float(np.sqrt(np.mean(distances ** 2))), "common": len(common),
            "coverage": len(common) / length}


def compare_structures(left, right, minimum_common=30):
    comparisons = [(compare_backbones(a, b, minimum_common), a["chain"], b["chain"])
                   for a in left for b in right]
    result, a, b = max(comparisons, key=lambda item: item[0]["score"])
    return dict(result, left_chain=a, right_chain=b)


def unique_backbones(chains):
    representatives = []
    for chain in chains:
        if not any(compare_backbones(chain, existing)["score"] > 0.995 for existing in representatives):
            representatives.append(chain)
    return representatives


def geometric_features(chains):
    descriptors = []
    for chain in chains:
        xyz = chain["xyz"]
        centred = xyz - xyz.mean(axis=0)
        eigenvalues = np.linalg.eigvalsh(centred.T @ centred / len(xyz))
        distances = squareform(pdist(xyz))
        separation = np.abs(chain["residues"][:, None] - chain["residues"][None, :])
        nonlocal_pairs = separation > 4
        descriptors.append({
            "ordered_residues": len(xyz),
            "radius_of_gyration": float(np.sqrt(eigenvalues.sum())),
            "span": float(distances.max()),
            "end_to_end": float(np.linalg.norm(xyz[0] - xyz[-1])),
            "contact_density": float(np.mean(distances[nonlocal_pairs] < 8)),
            "axis_ratio": float(eigenvalues[-1] / max(eigenvalues.sum(), 1e-12)),
            "planarity": float(eigenvalues[0] / max(eigenvalues.sum(), 1e-12)),
        })
    return {key: float(np.mean([row[key] for row in descriptors])) for key in descriptors[0]}
