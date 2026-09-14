# Certified numerical limits of the current spectator family

Independent exact replay now brackets the optimal periodic lower certificate
in the fixed current family to less than `2e-7` per site for both W=0 and W=1.
This resolves the numerical-limit question at that precision for these two
targets. It does not prove exact attainment of an optimum, a physical quantum
extension, or a solution to general quantum chemistry.

The authoritative new artifacts are in
`results/marginal_graded_hubbard8/spectator_hopping/<case>/thermal/final/`.
The root `thermal_limit_summary.json` records exact fractions and receipt paths;
`thermal_limit_provenance.json` records construction artifacts and source hashes.
Earlier directories remain historical evidence.

| Target | Accepted periodic lower | Accepted family ceiling | Exact interval width, decimal | Reduction from previous width |
|---|---:|---:|---:|---:|
| W=0 | −0.643058131811614 | −0.643057945542764 | 1.86268849433e−7 | 307.84× |
| W=1 | −0.660921390959772 | −0.660921298140620 | 9.28191515097e−8 | 4999.54× |

Each ceiling uses 77 positive physical source vectors, within the unchanged
85-source parser limit. Each exact mixture satisfies all 85 retained moment
equations and the fidelity constraints. The ceiling bounds all unrestricted
correction coefficients and reflected mean-correct range-two profiles in the
same fixed family: projector sources, ratio, ceilings and nine sparse shapes
are held fixed. A family ceiling limits the lower bounds this method can
produce; it is not an upper bound on the physical ground energy.

Both lower and upper family endpoints improved. The new
`interval_refinement.json` receipts check the fixed fields, audit the input
family receipts' hashes, and use rational arithmetic to verify nested intervals
and strictly smaller gaps relative to `completed_basis/`.

## Physical energy intervals

The targets remain half-filled open million-site U=4, t=1, V=1/2 chains, with
W multiplying the next-nearest density interaction in q=n_up+n_down−1.

| Target | Open-chain lower/site | Physical upper/site |
|---|---:|---:|
| W=0 | −0.643060631811614 | −0.610676347051188 |
| W=1 | −0.660925890959772 | −0.618424482369328 |

The physical upper endpoints are unchanged but were freshly replayed. The
lower replays cover all 4096 local Fock states through 94 symmetry blocks,
with maximum local PSD dimension 200. The physical upper replay recomputes
the physical source and compares exact 24-site contraction with a 160-bit
outward enclosure before the million-site transfer calculation. No production
verifier, physical gate or certificate version changed (ENERGY v13/FAMILY v9).

## Construction and failed proposals

Discovery replaced the nonsmooth minimum eigenvalue with the full-spectrum
soft minimum

\[
s_\tau=m-\tau\log\sum_{i=1}^{4096}e^{-(\lambda_i-m)/\tau},
\qquad m=\min_i\lambda_i.
\]

Every objective call checked the bound
`m − tau*log(4096) <= s_tau <= m`. Gibbs expectations provide gradients;
spectral covariance provides diagonal and then full Hessian preconditioning.
Two random gradient directions and selected curvature components were checked
against finite differences. These numerical checks guide proposals only.
Curvature truncation and floating point Hessian regularization are confined
to discovery and do not enter physical acceptance.

Five bounded trials per target were retained: initial annealing, polishing,
diagonal scaling, expansion of the 52 signed-charge search bounds from ±2 to
±8, and full-metric scaling. Every trial had at most 500 full-spectrum
evaluations. Several stages hit iteration limits; the initial runs reached the
evaluation budget. Expanding the artificial search box and using the full
metric materially improved W=1. No production coefficient limit was changed.
The final W=0 lower comes from `thermal/scaled/`; the W=1 lower comes from
`thermal/fullmetric/`. W=0's later continuous points did not yield a better
rounded lower certificate.

Optimizer success flags did not produce valid dual certificates. Gibbs atoms
at the final full-metric points still had maximum moment residuals about
`4.50e-4` for W=0 and `4.79e-5` for W=1. The eight dominant vectors were used
as untrusted physical candidates. Rounded and perturbed vectors, prior sources
and determinant anchors produced 4900 candidate columns, with zero pricing
rounds. The existing exact-weight construction succeeded directly for W=1.

W=0's selected 77 columns were inconsistent with the full 85-equation system.
The failed diagnostic remains in `thermal/final/`; the accepted proposal
instead comes from `thermal/affine_completion/`. The retained candidate ledger
and the previous accepted 85-source pool were reduced with exact integer
elimination to eight residual equations. Sixteen bounded residual LPs each
selected two old columns whose residual rank was one. A strict independent
column solver correctly refused those dependent selections.

The affine completion variant retained the free parameter and derived its
exact feasible interval from nonnegativity of both completion coefficients
and all original selected coefficients. Choosing an energy-minimizing endpoint
gave exact feasible solutions for all 16 attempts. The selected completion
uses old source indices 45 and 62; some original coefficients become exactly
zero, leaving 77 sources. It rechecks all original 85 integer equations and
all source, trace and rational-size limits before exporting a proposal.
The independent production family verifier then accepted it. No moment row
was dropped and no negative weight was clipped.

## Quantum consistency remains obstructed

Fresh independent CAR reconstructions confirm exact zero for all 14 spectator
moments, four spin moments and the unconditioned hopping moment in both new
mixtures. The averaged signed-charge law has an exact stationary classical
Markov extension. That construction does not extend the quantum density
matrix, its coherences or its fixed global particle number.

All four independent pair-transfer tests still fail in both new mixtures.
For `J(i,j)=d_i†d_j+d_j†d_i`, the tested operator is the translated difference
of `Y=J(i,j)−J(4−j,4−i)` on adjacent five-site windows.

| Pair | W=0 moment | W=1 moment |
|---|---:|---:|
| 0,1 | +0.00117286714593 | +0.000409688907511 |
| 0,2 | −0.000163280531762 | −0.000131657130934 |
| 0,3 | +0.000285884805336 | −0.000975262890186 |
| 1,2 | +0.0000145279070316 | +0.0000304883030521 |

The receipts store exact nonzero rational moments. Full-Fock Hermiticity,
PH/reflection/spin-flip symmetry, periodic cancellation, rank four modulo the
existing correction class, and a norm upper bound of four were rechecked.
These violations survive the stated averaging and exclude a stationary
quantum extension matching these particular averaged local density matrices.
They do not prove that every optimal witness violates these constraints or
that adding them will strictly improve the energy bound. Energy integration
of the pair-transfer directions remains the next unresolved experiment.

## Validation and remaining scope

There are 16 new accepting receipts: energy, family, interval refinement,
pair transfer, spectator overlap, spin overlap, coherent overlap and classical
Markov extension for each target. The summary audits these and the previous
81 receipts against current source hashes. The new interval comparison has
focused refusal tests for changed families, invalid gaps, rejected inputs and
non-nested intervals. Together with the unchanged fraction-free helper tests,
25 focused tests passed in 0.07 seconds. The full suite was not rerun; its
earlier 847 tests and 102 subtests remain historical validation only.

All work used the local CPU; observed run times are not controlled performance
benchmarks. No GPU or paid resource was used in this continuation. Construction
and verification are finite experiments, not a requested-accuracy scalability
theorem. No new ablation was performed at the selected lower coefficients.

W=+0.1 and W=−0.1 retain their previous `polished/` checkpoints; their numerical
limits were not tightened here. Transfer to generic molecular, long-range or
higher-dimensional systems, and general quantum representability, remain
unproved. The broader goal remains active.
