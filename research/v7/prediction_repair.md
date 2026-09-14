# v7 prediction request repair

The request gate accepts training design/output rows, observed history, proposed
inputs, and an explicit training residual enclosure. It never accepts future
outputs. For an input query outside the training row space, an exact nullspace
witness constructs two compatible affine ARX parameters whose predictions differ
by four times the requested tolerance. The gate abstains and returns that
witness. Row-space queries remain supported with the declared support bound
`eta ||w||_1`; rank deficiency alone is not treated as a refusal.

The result is conditional on the unconstrained affine raw ARX class selected by
the v6 fitting protocol. The residual enclosure measures fit/data uncertainty and
does not establish physical applicability.


## Actual measured-record replay

`python -m experiments.v7_measured_gate` reads only the TCLab step record.
The existing V6 fitter chooses ARX order 4. All 897 admissible lagged training
rows, including the available selection/calibration-era step observations,
enter the compatibility set. The proposed exact rational coefficient vector
has maximum residual approximately 0.6820845905 degrees C. This is a measured
residual enclosure, not an independently validated physical noise bound.

For current-Q1 coefficient index `p`, choose `v[p]=1`, `v[intercept]=-50`, and
zero elsewhere. Exact arithmetic verifies `Xv=0`. At a requested Q1=75,
`r^T v=25`. The constructed compatible parameters disagree by 0.4 degrees C,
exceeding twice the requested 0.1-degree tolerance. That one-step ambiguity
blocks the requested 16-step guarantee before any future response is read.
The artifact retains both parameter vectors, query, null vector and exact
training residual bounds so the witness can be independently replayed.

A query equal to a training row is in the row space despite nonunique parameters.
Its model-prediction interval has radius eta; that can itself exceed a requested
engineering tolerance. The gate distinguishes those two refusal reasons. An
arbitrary future query at the same input is not assumed informative merely
because one historical row is informative. Additional future process noise or
model misspecification requires a separate applicability/error allowance.
