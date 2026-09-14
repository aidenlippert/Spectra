# Compact diagonal overlap corrections

Eight symmetry-invariant five-site diagonal correction shapes break the
previous fixed-projector/profile ceiling with an accepted exact certificate. The periodic bound is
-0.6436859905704608 per site on the periodic U4,t1,V1/2 chain. The preceding
family's exactly proved upper limit is -0.643986740577249, so the accepted certificate
exceeds that limit by about0.0003007500. The open-chain lower is -0.6436884905704608 per site after the2.5/N
closing-bond correction. The exact replay took539.526 seconds; all482
regression tests passed in485.605 seconds. The historical physical upper
remains -0.6106763470511881. The resulting interval is about0.90% narrower
than the immediately preceding joint-projector interval.

## Which missing physical information is used

A translation-invariant physical state's five-site reductions of adjacent
six-site windows must agree. The previous exact local-mixture limit witness
does not satisfy this condition. Averaging it under signed spin exchange
and particle-hole symmetry preserves the Hamiltonian, the half projector,
and the charged projector family. Those source symmetries were checked
exactly rather than presumed from floating amplitudes.

Even after this averaging, the diagonal marginal difference at determinant
358 is about0.006381825035712739. Its reflection partner is613. Thus the
mismatch is not merely an arbitrary choice of a spin-asymmetric mixture.
The actual diagonal entries are rational and are stored in
`joint_projector/signed_density/symmetry_averaged_overlap_separator.json`.

The five-site occupation basis is grouped into spin-exchange/particle-hole
orbits. Each orbit paired with its spatial reflection yields a diagonal Y
that is invariant under those internal symmetries and odd under reflection.
Orbits that contain their own reflection are discarded. The eight shapes
with largest exact expectations in the symmetrized mixture were selected
for a bounded search. Their union has56 nonzero occupation entries.

## Why the correction preserves the target Hamiltonian

For any such five-site diagonal Y, define the six-site correction

    T_Y = Y(first five sites) - Y(last five sites).

Summing all periodic translates cancels every five-site term once with
positive and once with negative sign. No projector-window boundary term is
omitted. The usual open-chain correction still removes only the closing
Hamiltonian hopping and density bond.

Reflection reverses the two embeddings while sending Y to -Y, so T_Y is
reflection-even. It is diagonal and conserves particle and spin numbers.
Consequently the local all-Fock proof uses exactly the existing94 symmetry
blocks, maximum dimension200, and retains one rank-one projector update per
block. No PSD matrix cap was raised.

The v5 accepting path reconstructs T_Y from the supplied rational diagonal,
checks reflection oddness, and includes its exact orbit-norm-scaled diagonal
in the local positivity matrix. A new test independently checks the local
minimum on all4096 six-site determinants and periodic cancellation on all
65536 eight-site determinants. Refusal tests cover unsupported certificate
versions, malformed entries and corrections with the wrong reflection.

## Numerical construction and accepting artifacts

The search jointly adjusted six profile parameters, alpha, beta and eight
correction coefficients, with a hard250 matrix-evaluation budget. The U4,t1,
V1/2 run used70 evaluations and independently rebuilt the final CAR matrices.
The largest affine/fresh discrepancy was below1.3e-14. The rounded profile is

```
a=.278537 b=3.535971 p=.614080 q=1.155548
d=-.155883 e=.858821 alpha=.160840 beta=.121128
r=.5 ell=-3.0619549
```

The full56-entry correction is in
`results/marginal_graded_hubbard8/joint_projector/signed_density/symmetric_diagonals_8/profile_joint_r1_2_certificate.json`.
The accompanying `profile_joint_r1_2_certificate_replay.json` is the accepted exact
energy receipt. This is a compact
local consistency extension, not a general representability oracle.

## Transfer test

The same half/charged source vectors and the same eight correction shapes
were applied to U4,t1,V=-1/2. Only their scalar coefficients and local profile
were retuned; the numerical proposal is -0.5510766692912333 per periodic site.
This tests transfer when the density interaction changes sign within the
one-dimensional nearest-neighbor family. It does not cover generic molecular
integrals, longer-range terms or higher dimensions.

`discovery/joint_transfer_replay.py` independently replays the lower bound,
reconstructs the matched physical transfer state, checks a24-site exact
contraction inside its outward enclosure, and evaluates the million-site
upper. Its result is `symmetric_diagonals_8/transfer_V-1_2/matched_transfer_replay.json`.
The upper recipe was already known for this target; its physical expectation
is recomputed here rather than inferred from the positive-V result.

## A measured verification lead

The current exact PSD path uses fraction-free elimination. A separate
bounded probe on one actual132-dimensional tight joint Gram block proposed
an integer upper-triangular R and computed R^T A R exactly. Nonzero diagonal
entries prove R invertible. Exact strict diagonal dominance with positive
diagonal proves R^T A R positive definite and hence A positive definite.

The integer check took0.461 seconds; proposal0.022 seconds and Gram
preparation5.844 seconds in that local run. This is one block, not a controlled
whole-replay speedup. The witness and receipt are in
`joint_projector/congruence_probe/one_block_receipt.json` and its driver is
`discovery/gram_congruence_probe.py`. A separate standard-library replay rebuilt the actual Gram and used generic
integer matrix products; it also accepted, in3.893 seconds including
reconstruction. Full integration must preserve singular
PSD handling and all exact refusals. No GPU ran this probe.

General accuracy-versus-cost scaling and general representability remain
unproved. The old family cap does not cap this newly expanded family; a new
dual witness would need to cancel the added correction expectations too.

The transfer replay accepted in521.532 seconds. Its million-site interval
is [-0.5510791692912332, -0.5240746138143995] per site (displayed approximately),
closing 37.5384% of that target's preceding interval.
All source/input hashes in the receipt were checked after the run.
