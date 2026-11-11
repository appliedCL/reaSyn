# Development dataset: initial inventory

Inspected on 11 November 2026. This inventory covers the downloaded
`protein+ion_trajectories.zip` from de Bruyn et al. It distinguishes the contents
of this archive from the broader study described on the deposit page.

- Study: <https://doi.org/10.1021/acs.jcim.4c01172>
- Pinned deposit: <https://doi.org/10.5281/zenodo.13892825>
- Local archive: `raw/de_bruyn_2024/13892825/protein+ion_trajectories.zip`
- Size: 1,173,814,604 bytes.
- Published MD5, verified locally: `fa72beb77a32ae3325621d83f71c6422`.
- Local SHA-256: `cd39056ce4f045713ab4bbe49908960549979bd0e54514a9e3a6cffa51aedc8e`.

## Files and trajectory units

The archive contains two condition folders. Each holds `proteins+ions.pdb`,
`proteins+ions.gro`, `proteins+ions.tpr` and `proteins+ions.xtc`.

| Trajectory label | Archive folder | PDB/GRO/XTC atoms | Stored frames | Time span | Frame interval |
| --- | --- | ---: | ---: | --- | --- |
| AS | `DES-Amber simulations/AS/` | 2,384 | 60,001 | 0–600 ns | 10 ps |
| pAS | `DES-Amber simulations/pAS/` | 2,355 | 60,001 | 0–600 ns | 10 ps |

There are 120,002 stored frames in total. Coordinate uniqueness and statistical
independence have not been established. No a99SB-disp or pASH trajectory is present
in this archive, although those simulations are described in the broader deposit.
The other archives have not been inspected.

The paper reports analysing 100–600 ns. That inclusive interval contains 50,001
stored frames per trajectory, or 100,002 frames in total. Adopting this analysis
window remains a processing decision; the original archive is unchanged.

## Structural identity and condition metadata

| Field | AS | pAS | Evidence and variation |
| --- | --- | --- | --- |
| Protein residues | 1–140 | 1–140 | One C-alpha atom per consecutively numbered residue in each PDB |
| N-terminal cap | ACE at residue 0 | ACE at residue 0 | Present in both reference structures; retain separately from the 140-residue sequence |
| Residue 129 | SER | SEP, including a phosphorus atom | Phosphorylation state varies |
| Residue 50 | HIE | HIE | Same histidine naming in both structures; retain the original residue names |
| Force field | DES-Amber | DES-Amber | Folder labels and source description; no variation in the downloaded archive |
| Physical temperature | 300 K | 300 K | Source description; TPR parameters have not been decoded |
| Nominal salt | 150 mM NaCl | 150 mM NaCl | Principal-system conditions reported in the paper |
| Sodium/chloride atoms | 187 / 177 | 172 / 160 | Counted directly in the reference structures; ion counts alone do not establish a confounder |
| Water in selected files | Absent | Absent | Reference structures contain protein, acetyl cap and ions |
| Sampling method | REST2 | REST2 | Source describes the lowest-temperature state of 32 exchanging replicas |
| Simulation software header | GROMACS 2022.6, PLUMED 2.9.0 | Same | Version string in each TPR; full parameters remain unparsed |
| pH, run seeds and batch identifiers | Not established | Not established | Not recovered in this inspection |
| Original conformer population weights | Not established | Not established | No population-weight sidecar file is present in the selected archive |

Mapping HIE to H and SEP to S for sequence comparison gives the same 140-residue
sequence in both reference structures:

```text
MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA
```

Sequence SHA-256, using the ASCII sequence without a newline:
`9a06387610edd099466d614941c42a433f444feb6a8f25375bb0b6ad8ec9c6f4`.
This comparison preserves the distinction between sequence identity and chemical
modification. The PDB chain identifier is blank in both files.

## Hierarchy and identifiers

Preserve the hierarchy as study → deposited trajectory unit → stored frame.
Use `zenodo:13892825` as the pinned source accession. Within it, the archive member
path identifies the trajectory. A zero-based frame index, together with the
original time in ps, identifies one stored conformation within that trajectory.

For example, frame 10,000 of `DES-Amber simulations/AS/proteins+ions.xtc` is at
100,000 ps. Its reference must retain both the archive path and member path.
AS and pAS are condition/trajectory labels, not recovered seed or replica IDs.
The archive provides one trajectory per condition; independent replication has
not been established. Exchanging REST2 replicas must not be counted as independent
runs, and adjacent frames may be correlated.

## Checks performed

- The original archive's published MD5 matches the local file.
- All eight non-directory ZIP members were read successfully, including ZIP CRC
  checks, without extracting or duplicating the trajectories on disk.
- Every XTC frame header and compressed-payload boundary was traversed. Magic
  numbers, atom counts, finite times and box values, and strictly increasing frame
  times were checked. Both trajectories have a constant 10 ps interval.
- PDB and GRO atom counts and atom-name order match within each condition. Their
  reference coordinates are finite and agree after converting PDB Angstrom values
  to nm. XTC frame-header atom counts match those structures throughout.
- XTC coordinate scaling is 1,000 throughout. The file-format definition specifies
  coordinates in nm: <https://manual.gromacs.org/current/reference-manual/file-formats.html#xtc>.

Full trajectory coordinate decoding, TPR topology/parameter interpretation,
periodic-boundary handling and structural geometry checks remain outstanding.
Header consistency is not a complete coordinate or topology validation.

## Consequence for causal resampling

Phosphorylation state is a possible condition variable A. Force field, nominal
temperature and nominal salt do not vary across these two deposited trajectories.
No measured set of confounders C with sufficient independent observations has
been established. A field's availability or variation does not by itself justify
adjusting for it.

The current archive supports initial structural inspection and ordinary
resampling. It does not yet support estimating a non-constant causal correction.
With no justified varying confounders, the proposed stabilised rule
P(A) / P(A given C) reduces to unit weights. Equal trajectory representation then
gives equal selection probability to each condition; this does not establish a
causal effect. More frames from these same two trajectories do not supply more
independent runs.

Agree a causal graph, target distribution and suitable observations before
implementing an adjustment. Keep the reserved evaluation collection held out
while making those decisions. No causal weights or resampled dataset have been
created during this inspection.
