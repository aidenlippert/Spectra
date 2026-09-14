# Refined W=1 family ceiling and eighteen further consistency obstructions

The W=1 certificate-family gap has been reduced by **96.3062%**, from
0.005548129720482485 to **0.000204934492326946 per site**. The selected lower
energy certificate is unchanged. Exact physical replay accepts the new
118-source family ceiling. Independently, eighteen further three-spectator
consistency conditions fail in both selected local mixtures.

The authoritative W=1 directory is
`results/marginal_graded_hubbard8/two_spectator/W_plus_1/enriched_final/`.
W=0 remains `two_spectator/W_zero/final/`. The exact fractions, receipts,
construction diagnostics and source audits are in `two_spectator/refinement_summary.json`
and `two_spectator/refinement_provenance.json`.

| W=1 quantity, U=4, t=1, V=1/2 | Accepted value per site |
|---|---:|
| Periodic lower certificate | -0.6607111946158047 |
| Previous enlarged-family ceiling | -0.6551630648953223 |
| Feasible-box pricing ceiling, first pass | -0.6582903386310940 |
| Feasible-box pricing ceiling, refined pass | -0.6602065679378636 |
| Selected coherent-subspace ceiling | -0.6605062601234778 |
| Selected family gap | 0.000204934492326946 |
| Million-site open-chain lower | -0.6607156946158047 |
| Physical variational upper | -0.6184244823693281 |

Family ceilings bound the best lower certificate in the fixed family. They are
not physical ground-energy upper bounds. The remaining gap is nonzero and the
W=1 numerical limit remains unresolved. W=0's prior accepted family gap remains
0.000019009152447824. No physical energy lower bound or coupling-transfer result
changed in this refinement.

## Candidate-feasible pricing

The previous radius-clipping heuristic changed dual coefficients without
preserving the constraints for existing candidates. The replacement solves a
dual LP inside a radius-0.05 box while imposing every current candidate
inequality and the two fidelity-dual sign constraints. The box center moves to
the resulting solution. Only the pricing proposal uses this box; final family
verification has no extra box assumptions.

Before pricing, exact vector matching reduced the 14630-candidate ledger to
4211 distinct candidates while preserving all 118 physical sources of the
accepted mixture and all 4096 determinant anchors. Some accepted sources are
themselves anchors, so the counts overlap. This pruning preserves the existing
feasible mixture rather than assuming discarded candidates were unnecessary.
Ledger moment rows remain numerical proposals.

Forty rounds with up to two eigenvectors per physical block produced 4749
candidates and an accepted 117-source cap of -0.658290338631094. A second run
with eighty rounds and one eigenvector per block produced 5445 candidates and
the accepted cap -0.6602065679378636. Both passes price all 94 physical blocks.
The last reduced eigenvalues remained negative, approximately -0.04939 and
-0.01891. These runs do not establish pricing convergence. The numerical
candidate-feasibility residuals and every pricing round are retained in the
diagnostics. Observed run times were about 64 and 120 seconds on the local CPU;
these are not controlled performance benchmarks.

## Coherent subspaces and exact reconstruction

The first primal semidefinite problem allowed density matrices over the spans
of selected physical vectors, together with determinant anchors. It used seven
blocks of dimensions 1, 1, 39, 9, 15, 14 and 1. CLARABEL returned
`optimal_inaccurate`, with objective approximately -0.6604537282 and small
nonzero numerical residuals. This is not an accepted bound.

Four bounded conic rounds then enriched the spans using full-block reduced
eigenvectors. There were at most 32 active blocks, with a cap of 48 on each
compressed dimension; the largest reached dimension 45. Each solve was limited
to 30 solver seconds and 80 iterations. The block cap prevented some requested
enrichments, which are recorded. All four statuses were `optimal_inaccurate`.
The final numerical objective was approximately -0.6605066556, with equality
residual about 4.53e-8 and minimum numerical density eigenvalue about -2.25e-8.
The corresponding global reduced eigenvalue was still about -3.30. Thus neither
the conic status nor its dual certifies global optimality or exact feasibility.

Rounded physical eigenvectors and deterministic small perturbations were
added to the retained candidate ledger. The enriched attempt added 400 vectors,
giving 5845 candidates. An exact rational basis reconstruction produced 118
positive physical sources and one fidelity slack. Independent production
replay accepted all moment, normalization, fidelity and energy constraints.
The longest weight string is 2534 characters, below the existing 4096 limit.
The production family source cap stays 119; local energy verification is
unchanged. This accepted reconstruction gives the selected ceiling in the table.

The earlier, unenriched conic atoms exposed a separate conversion failure:
the numerical LP chose 118 columns for 119 exact rows. The original rectangular
rational solve correctly refused the inconsistent system. A bounded completion
helper computes an exact left-null relation, tests supplied physical columns
against its residual, and tries at most four nominees. A numerical QR step only
nominates independent rows; exact rational identities establish the relation
and the completed constraints. Negative weights remain forbidden.

The first nominee, a determinant anchor of mass approximately 2.8330e-12,
completed all 119 equations with nonnegative weights. This separate 118-source
physical mixture passed energy and family replay, giving a ceiling of
-0.6604538639170363. It is weaker than the selected enriched ceiling, but verifies
the exact completion method on the real failed conversion. The original failed
attempt is preserved in `W_plus_1/primal_atoms/` and its successful completion
in `W_plus_1/primal_completed/`. No residual was dropped or accepted by tolerance.

## Eighteen independent three-spectator obstructions

For a hopping pair i,j on five sites, use all three other sites k,l,m as
spectators and form `C = q_k^p q_l^q q_m^r B(i,j)`, with powers in {1,2}.
Particle-hole invariance requires `p+q+r+j-i` to be odd. There are forty such
primitive products. Four are reflection-invariant and vanish on
antisymmetrization; the other thirty-six give eighteen reflection pairs.
Set `Y = C - reflected(C)` and `T = Y_left - Y_right` on six sites.

The new probe reconstructs all eighteen CAR actions on all 4096 states,
checking Hermiticity, particle/spin sectors, reflection, particle-hole and spin
flip, plus periodic telescoping cancellation. A separate bit-swap implementation
checks every action. Each charge power has norm at most one and the disjoint
spin-summed hopping has norm at most two. The four-term telescope therefore has
norm at most eight.

The exact integer Gram matrix on two-bit off-diagonal entries with particle
number at most four has rank 47 for the preceding hopping directions and rank
65 after adding these eighteen. Fixed projectors vanish on that restriction;
the diagonal, spin and pair-transfer operators have no such entries. This proves
all eighteen new directions are independent modulo the current family.

All eighteen moments are nonzero in the selected W=0 and W=1 mixtures. Their
largest absolute moments are approximately 0.00011186450 and 0.00052225418,
respectively. Dividing by the norm bound gives normalized violations of at
least 0.00001398306 and 0.00006528177. Exact rational moments are in each
`three_spectator_overlap.json` receipt.

The thirty two-spectator, fourteen one-spectator, four pair-transfer, four spin
and coherent hopping moments still vanish in the selected W=1 mixture. Its
averaged charge law still has a classical Markov extension. The new nonzero
symmetry-invariant moments nevertheless forbid a stationary extension of that
particular averaged quantum mixture. They do not invalidate its family-ceiling
certificate, and no three-spectator energy improvement has yet been proved.

## Validation and remaining scope

Production sources and the previous 381-file provenance are unchanged. The
preceding full result of 945 tests and 102 subtests still applies to those
production sources. Seven new focused discovery tests passed separately in
2.58 seconds, covering all three-spectator actions and parity classes, exact
completion below floating-point precision, and infeasibility/negative-weight
refusals. This is not a claim that a new full 952-test suite was run.

New current receipts and their source hashes were audited, and the preceding
26 current receipts were re-audited without being represented as newly executed
replays. Historical source snapshots and earlier transfer evidence remain intact.
All launched numerical, exact-replay and test jobs are terminal. No GPU, paid
resources or additional agents were used.

The next mathematical work is to integrate or otherwise exploit the eighteen
new physical constraints, while continuing to address the remaining W=1 family
gap. Exact attainment, general quantum representability, generic molecular,
long-range or higher-dimensional transfer, and cost at requested accuracy
remain unproved. The goal remains active.
