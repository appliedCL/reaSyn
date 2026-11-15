# reaSyn

Resampling alpha-synuclein structural ensembles from the [Amyloid Atlas](https://people.mbi.ucla.edu/sawaya/amyloidatlas/) and PDB.

The workflow validates frozen inputs, measures backbone similarity, freezes dependency-aware temporal splits, runs six descriptive model baselines, exports capped resampling draws and reports uncertainty. The causal estimator is tested on synthetic confounding. Real-data causal correction remains unavailable because a defensible causal design and sufficient independent-sample metadata have not been established.

