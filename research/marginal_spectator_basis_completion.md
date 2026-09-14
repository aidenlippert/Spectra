# Exact completion of the deficient spectator-family bases

Both previously rejected numerical proposals now have nearby **exact 85-source
family certificates**. Completing the W=0 basis required two old physical
sources with weights of order `1e-12`; completing W=1 required one determinant
source with weight about `2.80e-10`. Every original moment equation, fidelity
inequality and positivity requirement was retained.

Independent production replay accepts both ceilings. This resolves the exact
basis-construction failure for these proposals. It does **not** resolve the
unrestricted relaxation's numerical optimum.

## Accepted results

The authoritative new artifacts are in
`results/marginal_graded_hubbard8/spectator_hopping/<case>/completed_basis/`.
`basis_completion_summary.json` at the spectator-hopping root contains the exact
fractions, provenance and receipt paths. Earlier `limit_refined/` and `polished/`
directories remain historical evidence.

| Target | Accepted periodic lower | New family ceiling | Remaining family gap | Reduction from preceding gap |
|---|---:|---:|---:|---:|
| W=0 | −0.643060225064887 | −0.643002883917473 | 0.0000573411474140 | 15.56% |
| W=1 | −0.661004875948134 | −0.660540822812267 | 0.000464053135867 | 9.20% |

The ceilings improve by about `1.05660254e-5` and `4.70179514e-5` per site.
W=0 now uses 85 mixture sources, compared with 77 in the preceding selected
witness; W=1 remains at 85. Both stay within the existing FAMILY v9 source cap.

The accepted lower certificates did not change. Their million-site open-chain
intervals remain:

| Target | Open lower/site | Physical upper/site |
|---|---:|---:|
| W=0 | −0.643062725064887 | −0.610676347051188 |
| W=1 | −0.661009375948134 | −0.618424482369328 |

These are the same half-filled U=4, t=1, V=1/2 targets. The lower and physical
upper were freshly replayed in each selected directory, including the full
4096-state lower check and the exact-24-site/160-bit upper enclosure comparison.
W=+0.1 and W=−0.1 retain their previous `polished/` checkpoints.

Each `ceiling_improvement.json` matches the fixed projector sources, ratio,
ceilings, sparse span, target, range-two profile and accepted periodic lower,
then proves a strictly smaller exact family ceiling. The comparison covers the
same unrestricted correction family. A family ceiling limits attainable lower
certificates; it is not a physical ground-energy upper.

## The repair

The saved finite candidate ledgers contain 8346 W=0 vectors and 8820 W=1
vectors. Their numerical LPs selected only 83 and 84 columns, respectively,
for a system with 85 retained equations. Previous exact solves found these
selected bases inconsistent. Row rescaling and adding fidelity-slack columns
had not produced a nonnegative exact solution.

The new approach keeps those selected columns and uses the 85 sources from a
previously accepted, strictly positive witness as a completion pool. It tests
every pair of pool columns for W=0 and every single pool column for W=1. No
spectral sampling or pricing is repeated.

Physical source columns are reconstructed independently from integer state
vectors. Multiplying a normalized source column by its integer norm removes
its normalization denominator. A common scaling of the two fidelity rows
clears the remaining fixed denominators, including the right-hand side. The
selected moment matrix, completion columns and target are therefore integers.

For selected columns \(A\), old completion columns \(B\), and target \(b\),
fraction-free elimination reduces the system to

\[
D x + Pz = r,\qquad Cz=s.
\]

There are only two residual equations for W=0 and one for W=1. Each possible
completion solves those residual equations exactly, then reconstructs the
selected coefficients. The algorithm checks nonnegativity of **all**
coefficients and checks the original 85 equations again, before exporting
normalized mixture weights. It also enforces the existing source count and
rational-weight length bounds. No negative weight is clipped and no equation
is omitted.

`fraction_free_completion.py` checks every integer division for an exact
remainder of zero. Work is bounded to 85 rows, at most 85 completion columns,
4000 enumerated completions, 4096-bit input integers and a 50000-bit pivot
budget. These are construction limits for this experiment, not a general
scalability theorem.

| Target | Selected columns | Completions tested | Nonnegative completions | Selected additions |
|---|---:|---:|---:|---|
| W=0 | 83 | 3570 | 14 | Old columns 56 and 81 |
| W=1 | 84 | 85 | 1 | Old column 39: determinant `21` |

The selected added normalized weights are approximately:

- W=0: `1.37469113e-12` and `3.06882464e-12`.
- W=1: `2.80444866e-10`.

The resulting exact ceilings are about `1.63e-11` and `3.79e-10` above the
previous numerical LP objective values. Thus very small physical components
repair the proposals, but the numerical values themselves were never accepted
as exact bounds. W=0's completion also enlarges the original finite candidate
pool with old source vectors; this is not a proof of optimality over the
original ledger.

The initial W=1 attempt stopped because the old witness included individual
determinants that were not represented by a single invariant-basis block. The
anchor variant evaluates these sources directly, including their spin-diagonal
moments. It then checks the entire old mixture's moments against the accepted
target. The failed log is retained as `completion_before_anchor.log`; the
corrected run and the independently accepted production replay are separate.

Observed preparation times were about 15 seconds per target. Integer reduction
took about 14 and 18 seconds, with final pivots of 4060 and 4284 bits. Total
construction times were about 31 and 34 seconds. The jobs ran concurrently;
these are observations, not controlled performance comparisons. Final weight
strings have maximum lengths 2429 and 2497, below the unchanged 4096-character
parser limit.

## Physical consistency remains incomplete

Fresh independent probes confirm exact zero moments for all 14 spectator
directions, the four spin directions and the unconditioned hopping direction.
The averaged charge laws still have the checked stationary classical Markov
extensions. Neither fact proves stationary quantum representability.

All four independent pair-transfer obstructions remain nonzero:

| Pair | New W=0 moment | New W=1 moment |
|---|---:|---:|
| (0,1) | +0.00120101580009 | +0.00100778585341 |
| (0,2) | −0.000206309477584 | +0.000198152522406 |
| (0,3) | +0.000236198578270 | −0.00112531139599 |
| (1,2) | −0.0000269587056373 | −0.000213831781884 |

The exact CAR, symmetry, rank-four independence, periodic cancellation and
norm-bound checks were rerun. These moments exclude a stationary quantum
extension matching the averaged six-site density matrix. They do not exclude
every inhomogeneous extension of one marginal. Pair-transfer corrections have
not yet been integrated into the energy certificate.

## Validation and remaining work

Nine focused tests pass for the new integer-elimination helper. They check
original-equation reconstruction, row swaps, a negative determinant,
inconsistent rectangular systems, explicit residual completion, and selected
refusal cases. Both exported witnesses then passed the unchanged production
physical verifier.

There are 16 new accepting receipts. The audit covers **81 current and
historical receipts with 2650 matching source-hash entries**. The production
ENERGY v13 and FAMILY v9 verifiers were not modified. The prior full regression
of 847 tests and 102 subtests was not rerun; the nine focused tests and fresh
physical replays are this turn's verification evidence.

No local calculation or GPU remains running. No GPU was used. The ledgers,
completion diagnostics, exact proposals and source hashes are retained for
continuation.

The remaining family gaps are still substantial, especially W=1. The numerical
optimum has not been proved. The successful basis repair can support future
candidate searches, and the four pair-transfer directions remain a concrete
route to stronger certificates. No new molecular, long-range or higher-
dimensional transfer was established. General representability and
requested-accuracy scalability remain unproved. The broader goal remains active.
