# Compact hopping consistency and remaining spin obstruction

For the open half-filled chain
`H = U sum doublons - t sum nearest hopping + V sum q_i q_(i+1) + W sum q_i q_(i+2)`,
with `q=n_up+n_down-1`, `U=4`, `t=1`, `V=1/2` and one million sites,
a compact off-diagonal constraint improves all four accepted lower bounds.
The physical upper bounds were freshly recomputed and are unchanged.

| W | Exact lower/site, decimal | Physical upper/site, decimal | Gain over signed-charge bound | Periodic lower-to-family-ceiling gap |
|---|---:|---:|---:|---:|
| 0 | -0.643181744195110 | -0.610676347051188 | 3.40995416e-07 | 1.95806101e-05 |
| 1/10 | -0.643828715474063 | -0.611451160583002 | 6.11559234e-07 | 1.2494495e-05 |
| -1/10 | -0.642802907793760 | -0.609901533519374 | 4.12909806e-07 | 1.3537428e-05 |
| 1 | -0.661623752383821 | -0.618424482369328 | 1.41963368e-05 | 3.1035482e-05 |

Decimals above display exact rational certificates. The family ceiling limits
attainable **lower certificates** in the stated family; it is not a physical
energy upper. The best current original-model interval is therefore
`[-0.6431817441951097, -0.6106763470511881]` per site. Complete rational values
and receipt paths are in [combined_summary.json](/Users/aidenlippert/Documents/Spectra/results/marginal_graded_hubbard8/hopping_telescope/combined_summary.json).

## What the new certificate adds

Let `B(i,j)=sum_spin(c_i† c_j+c_j† c_i)`. On five sites choose
`Y=B(0,3)-B(1,4)`; the six-site difference is
`T=Y_left-Y_right=B(0,3)-2 B(1,4)+B(2,5)`.
Every translated copy of Y cancels on the periodic chain. The certificate
adds `gamma T` to the local positivity inequality and then opens only the
physical Hamiltonian chain. There is no physical third-neighbor hopping
and no extra auxiliary opening penalty.

Production energy version v11 reconstructs the twelve CAR words, checks
Hermiticity, spin-number conservation and reflection invariance on all4096
states, and projects the operator into the existing physical reflection
blocks. All94 local blocks, maximum PSD dimension200, total dimension4096,
and151 total PSD checks are retained. The sparse diagonal cap stays64.
The scalar gamma is bounded exact nonzero input; older versions reject it.
The previous complete52 signed-charge directions remain present.

The numerical chart has67 variables: two hopping profiles, two projector
penalties, nine sparse coefficients,52 signed-charge coefficients, gamma,
and the local lower. The fixed old diagonal profiles are contained in the
complete charge span. Hard250 matrix evaluations apply to each search.
SLSQP reports convergence for W0 and W±.1; W1 reaches240 iterations.
Fresh physical reconstructions agree with affine search matrices within
3.56e-15, but these numerical checks do not replace exact PSD acceptance.
All four candidates passed standard-library-only exact lower and physical
upper replay, including an exact24-site trial-state contraction inside its
160-bit enclosure before the million-site recurrence.

Exact ablation receipts contain physical integer vectors with negative
Rayleigh residual after deleting only gamma T, and nonnegative residual
when it is restored. The negative residuals are about -.000666, -.000721,
-.000636 and -.005381 for W0,+.1,-.1,1. Thus this term is necessary for each
submitted coefficient set. Improvements over the preceding accepted bounds
also include further coefficient optimization: strict separation from the
**fully reoptimized** signed-charge-only family is not proved.

## Exact dual limits and numerical failures

Family version v7 adds the exact zero moment of T to every prior profile,
sparse-diagonal, quadratic, indicator and complete signed-charge condition.
It permits at most67 physical integer-vector mixture sources, all67 used
in each selected witness. Positive exact weights, trace one, both fidelity
inequalities and all physical moments are independently checked. The cap
covers unrestricted real correction coefficients with the fixed projector
sources, ratio, ceilings and nine sparse shapes, and all reflected
mean-correct profiles. It does not cover other supports or spin telescopes.

The initial W0 dual pool yielded a gap around4.15e-5. A bounded physical-state
pricing procedure used LP multipliers to select additional lowest-eigenvector
candidates across all94 blocks for60 rounds. Its exact accepted W0 witness
narrows the gap to1.96e-5. A clipped multiplier experiment was weaker; its
receipts and diagnostics are preserved separately. Clipped pricing stopping
with no new states does not certify an optimum of the unrestricted family.

For the other targets the first unrestricted pricing proposals failed exact
rectangular reconstruction, despite numerical LP success. These failures
remain in `unstabilized_diagonal_family_limit_diagnostic.json`; no new
proposal from those failed runs was accepted. Tightening LP tolerance from
1e-9 to1e-10 and retaining positive weights down to1e-14 instead of1e-10
produced67-column exact proposals. All three passed independent physical
replay and improved their family ceilings. Both numerical settings changed,
so the separate cause of recovery is not isolated. Selected candidates
range from6952 to7088 physical states; each run retains the60-round budget
and the13000-candidate hard cap. There is no conditioning, convergence-rate,
or large-scale performance claim.

These exact gaps remain around1e-5. The current family optimum is not sharply
resolved; optimizer success is not treated as its proof. The failed conic
cross-check from the preceding milestone was not rerun, and no GPU ran.

## Charge extension persists; spin extension fails

Independent exact probes confirm `<T>=0` in every new mixture. After explicit
PH/reflection averaging, each signed-charge law again has an exact stationary
order-five Markov extension: all243 prefix/suffix equations, stochastic rows,
stationarity and all729 reconstructed probabilities hold. This extends only
the averaged classical charge distribution, not the quantum density matrix,
spin law or a fixed global particle number.

Set `A(i,j)=4 S_i dot S_j` and
`Y_ij=A(i,j)-A(4-j,4-i)`, with local difference `Y_left-Y_right`.
All four tested spin telescopes still have nonzero exact moments:

| W | pair(0,1) | pair(0,2) | pair(0,3) | pair(1,2) |
|---|---:|---:|---:|---:|
| 0 | -0.006237601282 | -0.009899649131 | 0.0209310316 | 0.004284271487 |
| 1/10 | -0.01275315116 | -0.007188067556 | 0.01396815735 | 0.004655292189 |
| -1/10 | 0.001533630904 | -0.01303562407 | 0.02815606312 | 0.003823412589 |
| 1 | -0.03593976759 | -0.001642807368 | -0.02586575702 | 0.00940111335 |

The spin probe checks fresh CAR Hermiticity and PH/reflection/spin-flip
invariance on all4096 states, periodic coefficient cancellation, and an
independent two-site singlet/triplet block. Each pair operator has norm3;
triangle inequality gives a conservative local telescope norm upper bound12.
It does not claim the exact telescope norm. The nonzero moments survive
PH/reflection averaging, so these averaged quantum mixtures still cannot be
the six-site marginal of a stationary quantum chain. This does not rule out
an inhomogeneous extension of a single marginal.

The next concrete certificate attack is the four compact spin-dot telescopes,
with exact symmetry and positivity reconstruction before optimization and a
new matching dual family. None has yet been integrated into an energy proof.
Improving the existing family-limit numerical gap is also still open.
General quantum representability, arbitrary molecular/long-range/higher-
dimensional transfer, and requested-accuracy scalability remain unproved.
Finite chain successes are not a solution of general quantum chemistry.

## Validation and provenance

Focused tests:18 passed in25.82s, including independent bit-swap CAR actions,
all-block image reconstruction (no discarded leakage), periodic cancellation,
exact acceptance and negative-PSD refusal, old-version refusal, malformed
coefficients, prior dual rejection, new moment acceptance and source-cap gates.
Full repository regression:819 tests and102 subtests passed in491.19s.
One unrelated warning reports a returned tuple in the untouched
`tests/test_v4_calibration.py::test_world`; it is not a failed assertion.
The earlier596-test log covered only marginal tests; this run covered the
whole repository, so the totals are not a direct count of added tests.

All24 selected accepting receipts (six per target: energy, family, ablation,
hopping moment, spin moments and classical extension) have768 current source
and input hash entries, all verified. Changed pre-stage proof sources are
saved under `pre_hopping_sources/manifest.json`. Numerical candidates never
bypass the exact accepting verifiers. All research and validation jobs are
terminal; no GPU resource was started. The broad goal remains active.
