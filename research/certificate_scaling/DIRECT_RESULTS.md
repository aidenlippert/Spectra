# Direct certificate discovery: completed second campaign

12 September 2026. The strongest result is an **exact obstruction to one
plausible sparse-discovery family**, rather than a new scalable chemistry solver.
The campaign also completed 39 Lambda jobs, independently replayed all 39
resulting molecular intervals, and established a warm two-host research pool.

## A precise reason to stop tuning equal-weight pair squares

For square H4, consider all quadratic dictionary words from the existing
`dictionaries(8, 'quadratic')`, every singleton square and every equal-magnitude
signed two-word square, an arbitrary number multiplier of body order at most
two, and the full coefficient-l1 residual bound. An exact rational dual gives

```
maximum possible certified lower endpoint in this family <= -4.459007306545172 Ha
independently certified physical ground lower endpoint    >= -3.330388761351986 Ha
unavoidable error allowance                               >=  1.128618545193186 Ha
```

This is not an optimizer-convergence observation. Let y be the saved dual
functional. Exact checks give y(I)=1, y((Nhat-N)X)=0 for every allowed multiplier
basis element, y(B†B)>=0 for every declared square atom, and ||y||infinity<=1.
Consequently, for any certificate H=bI+SOS+(Nhat-N)X+R in this family,

`b - ||R||1 <= y(H)`.

The independent stdlib checker verifies 18,496 square constraints and 443
multiplier constraints, including generated degree-six words. The construction
also includes cross-charge pair atoms, a larger cone than the same-charge
implementation, so the obstruction applies to the latter as well. Missing dual
coordinates are explicitly zero. The physical comparison uses an already
verified lower certificate, not an FCI estimate.

The obstruction does **not** apply to arbitrary two-word PSD blocks, higher
degree words, larger multiplier families, different orbital/operator bases, or
general SOS. A same-charge two-word block has exact negative dual determinant
`-1/64`, exposing one direction this equal-weight cone misses. This is an
algebraic lead for better pricing, not proof that adding one atom closes the gap.

Reproduce without numerical packages:

```sh
python3 -S research/certificate_scaling/direct_dual_replay.py
```

The discovery script is `direct_dual_obstruction.py`; the proof and independent
receipt are in `results/certificate_scaling/direct_dual_obstruction/`. Discovery
currently constructs dense intermediate LP arrays before converting to sparse
matrices. That implementation is a bounded obstruction experiment, not a
memory-scaling achievement. The exact replay does not use that LP or its solver.

## Direct H-only discovery versus accuracy

The new LP chooses nonnegative sparse square weights jointly with b and the
number multiplier. It optimizes `b - ||R||1` and exports rational factors.
Only H, mode count, and particle count enter discovery. It never reads an
existing certificate's factors, supports, energy endpoints, or upper witness.
The selection rule ranks same-charge word pairs by their coefficient overlap
with H. This is a heuristic; zero overlap is not a safe exclusion theorem.

Thirty Lambda runs tested budgets 0, 128, 512, 2,048, and 8,192 with multiplier
body order one and two on square H4, rectangle H4, and the 12-mode H6 chain.
All pass exact lower replay. Increasing the budget beyond the generated pool
does not add new directions; those repeated settings are recorded, not counted
as different mathematical dictionaries.

Nine further runs compare H-ranked two-word PSD blocks with full quadratic
Gram matrices. The following are **complete independently replayed intervals**:

| Model | Fixed signed-pair LP, budget 8,192, body-2 ideal | Pair PSD, budget 1,024 | Full quadratic SDP |
|---|---:|---:|---:|
| H4 square, eight spin orbitals | 1.12863482 Ha | 1.12769051 Ha | 0.00503664 Ha |
| H4 rectangle, eight spin orbitals | 0.79705489 Ha | 0.46421533 Ha | **0.00085708 Ha** |
| H6 chain, twelve spin orbitals | 7.05718602 Ha | 27.74616871 Ha | 0.01302935 Ha |

The pair-PSD budget is smaller than the largest LP budget; this table is a
budget/quality comparison, not a claim that one cone uniformly dominates
another. Two pair-PSD solves report `optimal_inaccurate`; their entire exact
residual is charged. Only the full-quadratic rectangle baseline passes the
fixed 0.0015 Ha target. The compact direct methods do not pass it.

The full-quadratic baseline is itself H-only discovery, but uses the complete
quadratic Gram blocks and is not a new scaling result. Its complete interval
certificate sizes are 25,251 bytes (square), 19,694 bytes (rectangle), and
127,594 bytes (H6). They have different accuracy from the earlier higher-degree
certificates, so byte counts alone are not a compression victory.

All intervals concern the supplied rational finite-basis Hamiltonians. Basis,
continuum, and numerical integral-algorithm errors remain outside their scope.
Upper witnesses are separate validation controls: H4 uses 70 fixed-N entries;
H6 uses 200 nonzero entries from the existing FCI-derived control. They are
attached only after discovery. Exact sparse Hamiltonian actions recompute the
upper endpoint; this does not make witness generation scalable.

Reproduce all returned interval checks:

```sh
python3 -S research/certificate_scaling/direct_results_replay.py \
  --input results/lambda_runs/direct_discovery_pool/retrieved \
  --output results/certificate_scaling/direct_intervals
```

## Other parallel attacks and what remains unresolved

The locality-only SDP produced valid H4 and H6 lower certificates but loose
bounds. H6 required 1,770 blocks, 23,300 Gram scalar entries, and 26,600
coefficient rows; its interval is about 6.33537 Ha wide. A small largest block
does not guarantee a small or accurate total problem.

Direct rank-one factor optimization from H improved some crude baselines, but
failed to converge to useful molecular intervals. Its zero-SOS comparison did
not jointly optimize the ideal, so it is not evidence of an advantage over the
joint LP baseline. The molecular branch's original use of source supports was
removed; only corrected H-only runs may be used as evidence.

A constructive diagonally dominant quartic pair-Hamiltonian family has exact
lower=upper=0 certificates at 8, 16, and 32 modes. Actual CAR replay uses 9, 17,
and 33 rational factor rows. A half-filled determinant with no doubly occupied
pair supplies the zero-energy witness. The recognizer currently scans O(M²)
pair coordinates. This restricted frustration-free family demonstrates a
checkable sufficient condition, not applicability to generic Coulomb molecules.

The natural next attack is dual-guided unequal-weight factors and larger
operator supports, with explicit accounting for candidate generation and
remainder error. In particular, any pruning theorem must justify removing
operators whose Hamiltonian overlap is zero but whose dual square value is
negative. Increasing the current signed-pair budget cannot evade the exact
obstruction above. No efficient omitted-family improvement bound has been
proved for broad chemistry.

DSOS/SDSOS and SOS basis pursuit are established precedents, not project novelty:
[Ahmadi–Majumdar](https://doi.org/10.1137/18M118935X) and
[Ahmadi–Hall](https://arxiv.org/abs/1510.01597).
No comparison here establishes an advantage over Rubin–Low–DePrince.

## Compute and verification

Two A10 hosts are active, at $1.29/hour each, under the user's instruction to
keep them running between batches. The API credential is outside the repository;
the CLI replaces dashboard provisioning. See
`results/lambda_runs/direct_discovery_pool/README.md` for IDs, SSH commands,
source archives, lifecycle, and cost rate.

The LP batch ran four concurrent processes and finished in 15.05 seconds after
setup. The comparison batch ran three processes and finished in 26.55 seconds.
These jobs use CPU solvers; no GPU acceleration is claimed. Initial batches
ran at different times because the second machine was still being prepared.
Both environments are now available without another boot/install cycle.

All 119 downloaded output hashes were checked, and all 39 lower certificates
and their exact upper endpoints replayed locally under `python -S`. Six
discovery controls and six CLI controls pass. The dual obstruction has a
separate independent audit and stdlib replay. Early sparse-map and pair-map
prototypes were corrected before the frozen remote sources were used.

No production solver or verifier changed. The earlier interrupted chain work
and its outstanding full-suite verification remain unfinished and separate.
