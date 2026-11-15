# reaSyn

Resampling alpha-synuclein structural ensembles from the [Amyloid Atlas](https://people.mbi.ucla.edu/sawaya/amyloidatlas/) and PDB.

The workflow validates frozen inputs, measures backbone similarity, freezes dependency-aware temporal splits, runs six descriptive model baselines, exports capped resampling draws and reports uncertainty. The causal estimator is tested on synthetic confounding. Real-data causal correction remains unavailable because a defensible causal design and sufficient independent-sample metadata have not been established.

## Run

Use Python 3.14 and the existing downloaded data:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-analysis.txt
make test
make analyse
make notebook
```

Open `notebooks/analyse_fibril_structures.ipynb` for the results. Shared numerical code lives in `reasyn/` so notebooks, the command-line workflow and tests use the same implementation.

Configuration is in `config/fibril_analysis.json`. Results are written to `data/analysis/<run-id>/report.md` with metadata, similarities, frozen splits, model coefficients, predictions, resampling rows and diagnostics. Completed runs verify output checksums before reuse. The primary threshold and temporal boundary are exploratory defaults recorded before fitting. Do not tune them against reported test results and continue to call that test independent.

## Container

```sh
docker build -f Dockerfile.analysis -t reasyn-analysis .
docker run --rm --network none -v "$PWD/data:/workspace/data" -v "$PWD/data/raw:/workspace/data/raw:ro" reasyn-analysis
docker run --rm --network none reasyn-analysis python -m pytest -q
```

