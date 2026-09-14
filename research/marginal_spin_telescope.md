# Spin consistency strictly exceeds the previous hopping family

Four compact spin constraints now yield exact lower certificates above the
previous hopping family's certified ceilings, at all four tested targets.
This is strict separation from that entire fixed family, including all its
coefficient tuning. It does not invalidate its correctly scoped ceilings.

For the open half-filled chain
`H=U sum doublons-t sum nearest hopping+V sum q_i q_(i+1)+W sum q_i q_(i+2)`,
`q=n_up+n_down-1`, `U=4`, `t=1`, `V=1/2`, and one million sites:

| W | Lower/site | Physical upper/site | Gain over prior lower | Excess over prior periodic family ceiling | New periodic family gap |
|---|---:|---:|---:|---:|---:|
| 0 | -0.643155738782484 | -0.610676347051188 | 2.60054126e-05 | 6.42480255e-06 | 2.91969712e-05 |
| 1/10 | -0.643782808113971 | -0.611451160583002 | 4.59073601e-05 | 3.34128651e-05 | 1.382556e-05 |
| -1/10 | -0.642779738508018 | -0.609901533519374 | 2.31692857e-05 | 9.63185772e-06 | 2.51814285e-05 |
| 1 | -0.661259216539476 | -0.618424482369328 | 0.000364535844 | 0.000333500362 | 0.000163014237 |

The displayed decimals come from exact rational receipts. The original-model
interval is now `[-0.6431557387824844, -0.6106763470511881]` per site.
The family ceilings are limits on attainable **lower certificates**, not
physical ground-energy upper bounds. Full rational values and receipt paths
are in [combined_summary.json](/Users/aidenlippert/Documents/Spectra/results/marginal_graded_hubbard8/spin_telescope/combined_summary.json).
Final artifacts are in each target's `polished/` directory; parent directories
retain the initial accepted energy results and numerical diagnostics.

## The new physical constraint

Write `A(i,j)=4 S_i dot S_j`. The four canonical pairs are
`(0,1), (0,2), (0,3), (1,2)`. On five sites,
`Y_ij=A(i,j)-A(4-j,4-i)`, and on six sites `T_ij=Y_left-Y_right`.
Every translated Y cancels on the periodic chain. The local correction adds
no physical spin interaction and no extra opening penalty.

Production energy version v12 accepts one through four bounded exact nonzero
spin coefficients. It reconstructs the spin-dot CAR action, verifies
Hermiticity, spin-number conservation and reflection symmetry on all4096
states, and projects into the existing physical reflection blocks. Older
versions reject the new field. All94 local blocks, maximum PSD dimension200,
4096 total local dimensions,151 total PSD checks and the64-entry sparse
correction cap remain unchanged. Every prior charge and hopping correction
and the independent overlapping-projector ceilings are retained.

The exact rank receipt isolates double-spin-flip matrix entries in the
(1,1) two-electron sector. Their Gram matrix is
`[[32,0,0,-16],[0,32,0,0],[0,0,48,0],[-16,0,0,48]]`, with exact rank4.
Previous diagonal terms vanish on these off-diagonal entries; hopping changes
two bits rather than four. The actual fixed half/charged projectors occupy
six- or five/seven-electron sectors and vanish here. This proves four new
independent directions modulo the previous fixed-source certificate span,
including penalty variations. It does not cover arbitrary changed sources.

All four final energy candidates passed standard-library-only exact lower
and physical upper replay. Each upper was freshly evaluated, with an exact
24-site contraction inside its160-bit enclosure before the million-site
recurrence. Physical upper endpoints are unchanged. The same fixed projector
sources, ratio, ceilings, sparse span and trial-state recipe transferred
across W=0,+.1,-.1,1. This demonstrates this one-dimensional range-two
parameter transfer; it is not generic molecular or higher-dimensional transfer.

`family_escape.json` for each target freshly replays the previous hopping
mixture with current exact production code, matches both families' fixed
sources/ratio/ceilings/shape span and target, and compares its ceiling with
the accepted new periodic lower. All four exact differences are positive.
This proves the gain cannot be reproduced by tuning only the old family's
coefficients. Separate spin-ablation receipts give physical integer vectors
with negative Rayleigh residual after deleting only the spin correction at
the submitted coefficients, and nonnegative residual when it is restored.

## Numerical search and its unresolved limit

The enlarged chart has71 variables: two nearest-hopping profiles, two
projector penalties, nine sparse coefficients,52 complete signed-charge
coefficients, the existing hopping telescope, four spin telescopes and the
local lower. Initial searches hit240 iterations (241 matrix evaluations).
One bounded polishing pass was then performed for each target. W0 reports
success after95 evaluations and W-.1 after4; W+.1 and W1 still reach the
iteration limit. Each pass retains the hard250 matrix-evaluation cap, and
fresh physical matrices agree with affine search matrices within3.56e-15.
These numerical observations are not convergence or performance theorems.

Family version v8 adds all four exact zero spin moments to all prior
conditions. Its source cap is71, and every selected witness uses71 physical
integer vectors. Exact positive weights, trace one, both fidelity bounds,
all profile/charge/sparse/hopping/spin moments and physical energies are
checked independently. The cap covers unrestricted real spin coefficients
within this fixed-source family, not other supports or additional corrections.

Dual proposal searches use physical determinant anchors plus bounded
LP-directed eigenvector pricing:60 rounds per run, at most13000 candidates,
LP tolerance1e-10 and positive-support threshold1e-14. Initial candidate
counts for the selected W0,W-.1,W1 witnesses are7334,7540,7566. Their family
proposals were reused and freshly matched to the polished energy certificates,
since all coefficients covered by the dual may vary. For W+.1 the initial
70-column proposal failed exact rectangular consistency; no certificate was
accepted from it. Repeating the same bounded procedure at the polished seed
produced a71-column,7326-candidate proposal that passed exact replay. This
is numerical recovery, not a relaxation of any acceptance condition.

The enlarged family gaps remain about1.38e-5 to1.63e-4 per site. In particular,
the W1 limit is still loose. Neither exact numerical optimality nor a sharp
convergence/conditioning guarantee has been established. The old conic
failure was not rerun, and no GPU resource was started.

## What still fails to extend

Independent exact probes now find zero for the previously violated hopping
moment and all four spin moments. The PH/reflection-averaged signed-charge
law still has an exact stationary order-five classical Markov extension:
all243 prefix/suffix equalities, stochastic rows, stationarity and729
reconstructed probabilities are verified. This is only the averaged
classical charge law, not the quantum density matrix or fixed global N.

Nine tested charge-conditioned hopping telescopes nevertheless have nonzero
moments in every new witness. Let `p_k=q_k^2` project onto empty/doubly occupied
spectator site k, and `C(i,j,k)=p_k B(i,j)`, where B sums both spin hoppings,
`i<j` has odd distance and k differs from i,j. Use
`Y=C(i,j,k)-C(4-j,4-i,4-k)` and `T=Y_left-Y_right`.
The strongest tested moments, using a common conservative norm bound8, are:

| W | (i,j,k) | Exact moment, decimal | Lower bound on normalized violation |
|---|---|---:|---:|
| 0 | (0, 3, 1) | 0.00232599151056 | 0.000290748938821 |
| 1/10 | (0, 3, 1) | 0.00225723035774 | 0.000282153794717 |
| -1/10 | (0, 3, 1) | 0.00248308162887 | 0.000310385203609 |
| 1 | (1, 2, 0) | -0.0035739023646 | 0.000446737795575 |

The probe checks all4096 CAR actions for Hermiticity, spin-number conservation,
PH/reflection/spin-flip invariance and periodic cancellation. The norm upper
bound8 follows from a commuting spectator projector, hopping norm at most2,
and the triangle inequality; it is not the exact telescope norm. Nonzero
moments persist under PH/reflection averaging and forbid a stationary quantum
extension matching these averaged six-site density matrices. No claim is
made about extending a single marginal into an inhomogeneous state.

The next concrete attack is these nine compact conditioned-hopping directions,
starting with exact independence and symmetry checks, then production energy
and matching dual integration before optimization. No conditioned-hopping
energy certificate exists yet. Tightening the current spin-family limit is
also open. General quantum representability, arbitrary molecular/long-range/
higher-dimensional transfer and requested-accuracy scalability remain
unproved. The broad goal remains active; these finite results do not solve
general quantum chemistry.

## Validation and provenance

New spin tests independently construct local spin actions across all4096
states, verify projected images have no discarded components, check periodic
cancellation, exact safe acceptance and unsafe PSD refusal, older-version
refusal, malformed coefficients, old-witness rejection, new-moment acceptance,
and source-cap/mode gates. The focused spin+hopping run reports27 passes in
46.14s, including14 new spin tests. The earlier spin-only run reports12 passes
in16.45s, before adding the two family tests.

The full repository command completed successfully: pytest reports
`824 passed, 1 warning, 93 subtests passed in 548.25s`. A subsequent read-only
collection check lists833 test items; the run's summary counts are preserved
verbatim rather than converted into another tally. The warning is from the
untouched `tests/test_v4_calibration.py::test_world` returning a tuple.

All32 selected accepting receipts (eight per target) plus the rank receipt
have1107 current source/input hash entries, all checked. Changed pre-stage
sources are retained under `pre_spin_sources/manifest.json`, and selected
numerical seed/source provenance is recorded separately. Initial lower
certificates and failed dual diagnostics remain distinguishable from the
final polished artifacts. No local research or validation process remains
running, and no GPU was used.
