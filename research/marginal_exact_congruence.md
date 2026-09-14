# Exact congruence witnesses for the current Hubbard certificate

This records the earlier fixed support. Newer accepted support replacements
and transfer results are in `research/marginal_shape_selection.md`.

The nine-correction U4,t1,V1/2 certificate now has two accepting CPU replays.
The original fraction-free elimination accepted in401.725 seconds. A fresh
standard-library replay using optional integer congruence witnesses accepted
in15.498 seconds, reconstructing all151 blocks and reproducing the same
energy bound and every block's dimension, rank and nullity.

These are sequential local observations on the same certificate, with
different replay drivers. Their elapsed ratio is approximately25.92; this
is not a repeated controlled benchmark or an accuracy-versus-cost result.
Witness proposal took7.960 seconds separately. No GPU performed exact proof.

## What is proved

For each freshly reconstructed symmetric integer matrix A, the optional
witness supplies an upper triangular integer matrix R. Every diagonal entry
must be nonzero, so R is invertible. The verifier computes C = R^T A R using
integer arithmetic, then requires, for every row i,

    C_ii > sum_{j != i} |C_ij|.

Symmetry and strict positive diagonal dominance imply C is positive definite.
Congruence by invertible R then proves A positive definite. No approximate
eigenvalue, Cholesky result or input digest constitutes acceptance. The digest
only selects a candidate factor for the freshly reconstructed A; its exact
congruence is still checked.

Integer factors are bounded by10^30 in magnitude, dimensions by300 and
the supplied collection by256 factors. A malformed, singular, wrong-sized,
unused or insufficient factor is refused. Exact negative or zero dominance
margins refuse. Witnesses apply to positive definite blocks only.

Blocks without a witness retain the original fraction-free elimination,
including negative pivots, null-pivot coupling refusals and every exact
division check. The positive initial-divisor input gate also applies when
a witness is used. On the witness route no elimination is performed: the
separate congruence theorem supplies the proof directly.

## Cost and compactness

The file has73 unique factors, used by79 blocks;72 blocks used the original
fallback. Repeated identical matrices reuse a factor. The largest factor
has dimension200. The auxiliary JSON is4,082,986 bytes and contains647,212
integer entries, including triangular zeros. This is additional speed data;
the energy certificate remains independently replayable without it.

Both integer products take cubic arithmetic work in matrix dimension.
If entries of A have b bits and entries of R have r bits, intermediate
entry sizes grow approximately as b+2r+2log2(n), rather than as determinants
of progressively larger submatrices. This explains the approach's appeal
on these measured blocks, but does not establish an end-to-end complexity
guarantee. Near-singular matrices may need factors beyond the fixed bounds,
or may require the original fallback. No new physical constraint follows
from this numerical preconditioning technique.

## Artifacts and validation

Production entry: `experiments/marginal_projector_extendibility.py` accepts
an optional `psd_witnesses` keyword for v2/v4/v5 energy replay. The checker
is `experiments/marginal_congruence_psd.py`. Discovery and standard-library
replay modes are separate in
`results/marginal_graded_hubbard8/discovery/joint_congruence_replay.py`.
Proposal mode has an explicitly nonaccepting matrix-collection hook. Replay
mode never installs that hook and reconstructs the original physical proof.

Receipts and witnesses are under
`results/marginal_graded_hubbard8/joint_projector/signed_density/symmetric_diagonals_9/extra_31/`:

* `congruence_witnesses.json`: optional auxiliary proposal.
* `congruence_witnesses_replay.json`: accepted all-block exact receipt.
* `congruence_summary.json`: sizes, measured times and source hash checks.
* `profile_joint_r1_2_certificate_replay.json`: original exact elimination.

The original executed source versions are preserved with hashes under
`joint_projector/pre_congruence_sources/`. Twenty-two focused new tests pass,
including a non-diagonally-dominant SPD example, indefinite and singular
witness refusals, malformed inputs, stale factors, exact-divisibility
fallback, and a complete mixed witness/singular energy replay. Full regression
status is tracked in `results/marginal_final_validation.json`.

The certified million-site lower is
`-2011436115385613/3125000000000000 = -0.6436595569233962` per site.
The fixed nine-shape relaxation is now independently capped to within
6.58e-8 per periodic site; see `research/marginal_nine_shape_limit.md`.
General representability, generic molecular transfer and requested-accuracy
scalability remain open.
