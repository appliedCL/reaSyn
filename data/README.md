# Source datasets

Two alpha-synuclein studies are selected for initial inspection. Original files are
stored under `data/raw/` and excluded from Git. `sources.json` records the pinned
deposits, original filenames, expected sizes, published MD5 checksums and download
URLs. Source metadata are preserved as `record.json` beside each deposit's files.

The publication plan is to release processing methods, source references,
identifiers, resampling weights and results, while linking to the original
coordinate deposits.

## Development source: de Bruyn et al. (2024)

- Dataset: <https://doi.org/10.5281/zenodo.13892825>
- Paper: <https://doi.org/10.1021/acs.jcim.4c01172>
- Local directory: `raw/de_bruyn_2024/13892825/`
- Selected files: `README.md` and `protein+ion_trajectories.zip`.
- Selected size: 1,173,816,654 bytes, including the README.

The deposit describes full-length, N-terminally acetylated monomers, with
unphosphorylated and S129-phosphorylated forms. It includes DES-Amber and a99SB-disp
simulations, with a singly charged phosphorylated form listed only for a99SB-disp.
The deposited protein-and-ion trajectories represent the 300 K state of REST2
simulations. The paper specifies 150 mM NaCl for the principal systems.

The selected archive contains configurations and non-parameterised topologies
according to the deposit. The solvent-containing and additional-simulation
archives are not selected. Their details remain in `sources.json`.

## Evaluation candidate: Ollila and Malm collection

- Deposits: <https://doi.org/10.5281/zenodo.17607490> and
  <https://doi.org/10.5281/zenodo.17582035>.
- Peer-reviewed paper: <https://doi.org/10.1038/s42004-026-02155-8>.
- Local directories: `raw/ollila_malm/17607490/` and
  `raw/ollila_malm/17582035/`.
- Selected files: all 25 XTC trajectories and 25 accompanying TPR files.
- Selected size: 50,088,375,432 bytes.

These two deposits form one study collection. Filenames identify starting
structures 01 to 05, five force fields and 2,000 ns trajectories. A filename does
not establish statistical independence. The deposits describe protein-and-ion
trajectories without water. Compatibility between each trajectory's atom subset
and the accompanying TPR remains to be checked.

## Download verification and scientific status

Downloads use a `.part` suffix until their byte counts and published MD5 checksums
have been verified. Verified files retain their original names. A local
`download_receipt.json` records verification results and SHA-256 checksums.
The ZIP archive is kept unextracted to conserve disk space.

Downloaded bytes and checksums establish file integrity. Sequence, modifications,
solution conditions, independent runs, frame correlation and shared experimental
or fitting evidence still require assessment before fixing evaluation subsets or
estimating causal weights. REST2 replicas are coupled and must not be counted as
independent simulation runs.

The previous 9A1A versus PED00545 comparison remains development evidence. No
causal correction, resampling or model training has been performed here.
