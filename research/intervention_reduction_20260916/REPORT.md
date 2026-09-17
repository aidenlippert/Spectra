# Certified molecular dynamics and finite inverse control

This pass added an exactly checked forward-and-control capability on existing
molecular Hamiltonians. It did **not** solve general many-body computation or
establish a field-wide breakthrough. The useful new result is a reduced model
whose driven trajectory, observable changes, and a selected control objective
are certified in the original rational Hamiltonian.

The main unresolved dependency remains discovery and preparation: the successful
cases used all 400 H6 or 441 water configurations in their magnetic sectors.
The query checker subsequently discards those configurations and uses small
exact matrices. That is a measured reuse result, not cheap independent discovery.
The H8 dynamical test still fails. All previous ground-energy milestones remain
unchanged; this pass does not replace or invalidate them.

## What was accepted

Each case has its original interacting molecular Hamiltonian and two
noncommuting one-body controls: a spin-summed occupation difference and a
spin-summed hopping between two declared spatial orbitals. Each control is
bounded by 0.01 Ha. The horizon is 10 atomic time units, roughly 0.24 fs.

The initial state is the exactly specified normalized first integer basis
column. It is not asserted to be the exact molecular ground state, a pure
singlet, or an experimentally prepared state. The occupation is that of a
specified basis orbital, not an inferred atomic charge or measured chemical
yield. The pulses are Hamiltonian diagnostics, with no laboratory field mapping.

The response target was normalized state-vector error <=0.005 throughout the
trajectory. This implies <=0.01 absolute error for a spatial occupation in
[0,2]. The following bounds include every polynomial residual, omitted action,
rounding term, segment jump, normalization correction, and a pointwise control
perturbation allowance of 1/400000 Ha per channel, subject to the original box.

| Case and pulse | Directions | Uniform state-error bound | Target met |
|---|---:|---:|---|
| H6, two-stage switch | 16 | 0.004429554286 | yes |
| H6, interior constant control | 16 | 0.001935148156 | yes |
| H6, three-stage switch | 16 | 0.004051252120 | yes |
| Water, two-stage switch | 40 | 0.002234682018 | yes |
| Water, interior constant control | 40 | 0.002798269464 | yes |
| Water, three-stage switch | 40 | 0.004604193613 | yes |
| H8, best tested switch model, nominal controls | 32 | endpoint bound 0.578344537739 | no |

The three test pulses were not training trajectories. H6 uses one unchanged
rank-16 model for all three. Water needed adaptive residual enrichment from
rank 16 to 40. Its rank-32 model passed two pulses and missed the third. This is
limited adaptive transfer, not success of the original untouched frozen rule
or a statistical reliability study.

H6's three robust occupation-change intervals are, respectively,
[-0.0254858,-0.00776761], [0.00534122,0.0130818], and
[0.0103223,0.0265273]. Thus both signs of a response are established. The final
water intervals contain zero; those tests establish bounded predictions, not
a decisive signed population change.

The control-neighborhood result covers infinitely many measurable perturbations
around each nominal pulse. It compares the perturbed true evolution to the
nominal reduced prediction. It does not establish the 0.005 target over the
entire original control box, all preparations, or longer times.

## A constructed policy with a proved effect

Before the inverse search, the target was an increase of at least 0.025 in
H6 spatial orbital 3, with the same initial state, horizon, and control limits.
A bounded search evaluated 729 three-stage schedules in the existing reduced
model. It selected this sequence, where D is the occupation-difference control
and W the hopping control:

| Time interval (atomic units) | D amplitude (Ha) | W amplitude (Ha) |
|---|---:|---:|
| 0 to 3 | -0.01 | +0.01 |
| 3 to 7 | +0.01 | +0.01 |
| 7 to 10 | +0.01 | -0.01 |

The initial rank-16 check proved the functional target but missed the stricter
state-error tolerance. Adding eight residual directions, without changing the
pulse or initial state, produced:

- State-error bound including the standard perturbation allowance:
  **0.003770050390**, uniformly over the trajectory.
- Certified final occupation increase: **[0.068383495449, 0.083463697009]**.
- A subsequent larger perturbation allowance, up to **0.0005 Ha per channel**
  while staying inside the original amplitude box, still certifies an increase
  of **[0.028583495449, 0.123263697009]**. This is 5% of the declared maximum
  control amplitude and still clears the requested 0.025 effect.

The larger perturbation allowance has a state-error bound of 0.023670050390, so
it does not meet the separate 0.005 state-accuracy target. The functional claim
and the state-accuracy claim are deliberately distinct. Global optimality of
the selected pulse, physical preparation, and synthesis are not proved.

The machine-readable result is
[CONTROL_POLICY.json](/Users/aidenlippert/Documents/Spectra/results/intervention_reduction_20260916/CONTROL_POLICY.json).

## What makes the proof reusable

For state columns V and controls H_l, the constructor computes the exact joint
action Gram matrix of

\[
Z=(V,H_0V,H_1V,H_2V),\qquad J=Z^\dagger Z.
\]

Every trajectory residual is a small coefficient vector acting through Z.
Its integrated squared norm is therefore an exact quadratic expression in J.
The checker retains all cross terms before taking a norm. After construction,
the query path requires no configuration records, Hamiltonian-sector matrix,
numerical eigensolver, or floating-point integration.

H6 rank 16 uses a 64-by-64 joint matrix; water rank 40 uses 160-by-160. An
observable Gram matrix and a metric are also retained. Matrix size alone is
not total storage: the entries are arbitrary-precision integers/rationals.
The matrix is not a self-authenticating proof of its own many-body definition;
the original CAR reconstruction is performed and charged once per bundle.

The exact trajectory theorem, uniform-time normalization enclosure, robust
control extension, kernel construction, and their hypotheses are derived in
[MATHEMATICS.md](/Users/aidenlippert/Documents/Spectra/research/intervention_reduction_20260916/MATHEMATICS.md).
The underlying Galerkin, Duhamel, POD, Krylov, and residual-certification ideas
are established methods. No field novelty is inferred from this implementation.
Primary-source context and mathematical review are in
[SOURCES.md](/Users/aidenlippert/Documents/Spectra/research/intervention_reduction_20260916/SOURCES.md).

## Measured cost and dependencies

The matched H6 comparison ran the same three proposals through direct and
compiled exact checkers, including the same robust and uniform-time claims.
Every mathematical receipt field matched exactly.

| Matched exact verification work | Seconds |
|---|---:|
| Direct, three queries including their preparation | 6.2650 |
| Compiled preparation once | 0.7459 |
| Three compiled queries | 1.4900 |
| Compiled complete batch | 2.2359 |

This single-host measurement is a **2.80x improvement in this verification
batch**. It is not an end-to-end solver speedup or a comparison with every
optimized conventional solver. Timing variability was not statistically sampled.

The final water batch took 17.37 seconds, including 2.40 seconds of exact
preparation and 4.67-5.15 seconds per query. Stable Taylor propagation had
required 280 short segments on the rank-32 switch case. Exact checking of a
Chebyshev polynomial proposal reduced that same query from about 15.32 to
2.80 seconds at the target, while retaining its full residual test. This is a
proposer/checker improvement; it did not change the physical representation.

Discovery is a separate and unresolved cost. The snapshot constructor uses
MPS-guided configuration selection and selected-space eigensystems. Subsequent
enrichment applies sparse CAR actions and adds residual directions without an
enlarged dense diagonalization, but materializes every reached label:

| Case | Initially selected labels | Final labels reached/retained | Full magnetic-sector dimension |
|---|---:|---:|---:|
| H6 | 400 | 400 | 400 |
| Water | 417 | 441 | 441 |
| H8 | 512 | 4,900 | 4,900 |

The supplied MPS and integrals were inherited. Their archived process times,
including stage startup, are listed rather than silently assigned zero cost:

| Input | Integral generation (s) | MPS state stage (s) |
|---|---:|---:|
| H6 asymmetric | 1.0202 | 14.8631 |
| H8 cold | 1.3224 | 234.6881 |
| Water asymmetric | 1.1612 | 23.5170 |

Adding historical component times to current times would not constitute a newly
timed cold pipeline. Older ground-energy optimization and verification were not
inputs to this experiment; their project R&D costs remain in the original ledgers.

Before final documentation and focused-test reruns, the campaign ledger contains
**65 instrumented processes**, **372.551 process-wall seconds**, **360.689 CPU
seconds**, and **221.512 MB maximum single-child RSS**. It includes unsuccessful
accuracy attempts, repeated checks, and reference calculations. Eighteen valid
but too-wide certificate receipts missed their numerical target. Process exit
success is not counted as numerical success. These totals exclude reasoning,
editing, browsing, documentation, the audit, unmetered focused tests, and prior
R&D; those are not zero-cost claims. See
[ACCOUNTING.json](/Users/aidenlippert/Documents/Spectra/results/intervention_reduction_20260916/ACCOUNTING.json).

## Independent checks and retained failures

An independent numerical reference constructs the full magnetic-sector matrix
using a separate literal occupied-orbital ladder implementation, then propagates
with sparse matrix exponentials. H6, water, and H8 predictions stayed inside
their accepted bounds. The numerical reference itself is not an exact proof.
It also remains fast: after preparation, the H6 switch propagation took about
0.021 seconds and H8 about 0.375 seconds. No competitive advantage over those
small-system propagators is established by this pass.

Focused tests exercise nonorthogonal coordinates, omitted configurations,
mixed control terms, large removable phases, exact two-level propagation,
segment jumps, malformed inputs, wrong bindings, intermediate norm zeros,
local CAR control norms, and direct/compiled equivalence. Independent inspection
reviewed the mathematical bounds and corrected stale findings in the source note.
The final clean-process run passed 18 accepting-boundary tests and four
numerical-proposer tests. Their logs are retained with the results.

The low-energy and response-weighted bases failed the original all-waveform
uniform envelope, including a rank-64 H6 attempt. H8's trajectory bounds improved
from 1.28627 at rank 16 to 0.68339 at rank 24 and 0.57834 at rank 32, far short
of 0.005. Independent propagation found state error about 0.13285 at rank 32:
both representation error and conservatism of the bound are present. These
bounded searches are not impossibility proofs for their attainable families.

The water rank-32 fixed-step Taylor candidate became unstable; its enormous
accepted error bound is retained. Stable propagation fixed that numerical issue.
The second adaptation, residual enrichment, addressed the remaining model error.
Unused candidate proposals, including two unverified inverse-search runners-up,
are not counted as successes. Early snapshot discovery metadata mislabeled POD
column indices as selected eigenvector indices; the writer is corrected and
the old receipts remain preserved. The first column is an eigenvector proposal;
the other snapshot columns are POD directions.

## What this changes, and what it leaves open

Spectra can now execute this restricted loop on actual molecular coefficients:
construct a reduced dynamical model, propose a control, and prove its effect
with an explicit robustness domain. This extends the work beyond a static
ground-energy interval. It supplies neither laboratory-ready controls nor an
inverse molecular architect.

The exact fixed-linear-subspace obstruction in MATHEMATICS.md rules out one
overly strong demand: a nontrivial proper state subspace invariant under every
number-conserving one-body spin-orbital control. It is a standard irreducibility
fact, not a new hardness theorem or an obstruction for the two tested controls.
Adaptive, nonlinear, tensor-network, approximate, and observable-specific
representations are not excluded.

The strongest unresolved scientific task is constructing the necessary
intervention response cheaply without first reaching the full sector. H8 makes
that gap concrete. Finite-temperature behavior, longer dynamics, reliable
physical models, molecular inverse design, preparation policies, and experimental
advantage remain unestablished. Existing SOS energy certificates offer possible
energy-constrained residual bounds, derived in SOURCES.md, but that route was
not demonstrated by the accepted molecular trajectory results here.

Reproduction commands and exact source dependencies are in
[RUNBOOK.md](/Users/aidenlippert/Documents/Spectra/research/intervention_reduction_20260916/RUNBOOK.md).
