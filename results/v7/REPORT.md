# V7: prediction validity repaired; method acquisition remains unproved

The history-state candidate is retired for the V6 workload. V1–V6 remain
preserved. The new prediction gate is operational on the actual measured-data
model class. Many-body calculations now have an exact residual checker and
costed conventional baselines. No two-generation method-learning result is
claimed, and no reserved evaluation outcomes have been opened.

## Request-specific prediction repair

The existing V6 fitter selects an order-4 affine ARX model from the TCLab step
record. All 897 lagged step-record rows enter an exact compatibility check.
The gate reads no future response measurements.

With constant training Q1=50, changing the current-input coefficient by h and
the intercept by -50h preserves every training prediction, including transient
history terms. At a proposed Q1=75, the model predictions change by 25h.
The saved witness consists of two exact rational models within the same
training residual enclosure (approximately 0.6820845905 °C) that disagree by
0.4 °C on the requested one-step prediction. This exceeds twice the requested
0.1 °C tolerance and blocks the requested 16-step guarantee.

The admitted class is unconstrained affine ARX with a free intercept. This is
not a claim that constant inputs always prevent identification. Query-specific
row-space certificates remain available: if `r = X^T w`, predictions over
`||X theta-y||∞ <= eta` lie within `eta ||w||1` of `w^T y`. Such a bound can
itself be too wide for a requested guarantee. The measured residual enclosure
is not independent evidence that this class describes hardware or that future
process noise has the same bound.

Evidence: [replay and compatible models](prediction_requests.json),
[derivation and scope](../../research/v7/prediction_repair.md),
[request API](../../experiments/v7_prediction_requests.py).

## Many-body mathematics and independent checking

For a finite supplied Pauli Hamiltonian with nonnegative uniform local
 depolarization, the dual semigroup is positive and unital. For a piecewise
rational polynomial P approximating `exp(t G) O`, contractivity and variation
of constants give

```
||exp(T G) O - P(T)||∞
 <= ||O-P(0)||∞ + sum(jump norms) + integral_0^T ||P'(t)-G P(t)||∞ dt.
```

The checker recomputes every residual coefficient and verifies norm-partition
coverage, pairwise anticommutation, exact outward square-root bounds and the
requested horizon. It supports coefficient-wise integration and an exact
Bernstein-basis envelope. Numerical proposals export a rational polynomial;
that polynomial's residual includes its numerical approximation error.
This is established residual mathematics implemented as a checkable procedure,
not a new universal many-body theorem.

Unlike V3's strong-damping theorem, this checker is valid at the prospectively
fixed gamma values 0, 1/5 and 2. A valid formula does not imply every calculation
fits the budget. All development arms use the same 1/1000 error tolerance,
T=1/5 or 1/2, declared generators, observables and sparse-support limits.
Independent dense two- and three-qubit calculations sanity-check the proof
implementation, including nonnormal damped dynamics. Larger cases are checked
through explicit sparse residuals, without a hidden full matrix propagation.

## Development calculations and rejected mechanisms

| Calculation or candidate | Certified development cases | Finding |
|---|---:|---|
| Coarse-grid full Taylor | 24/36 | Preliminary conventional comparator |
| BFS projected Taylor | 27/36 | Best feasibility count in this run |
| Residual-driven sparse expansion | 24/36 | No acquired procedure |
| Arnoldi projected polynomial | 24/36 | No acquired procedure |
| Strengthened adaptive Taylor | 26/36 | Standard recurrence reuse and adaptive order incorporated into baseline |
| Magnitude-based coefficient pruning | 26/36 | Conventional comparator |
| Work-weighted coefficient pruning | 26/36 | Slower on all 26 shared successes; rejected |

The first 144 calculation rows include construction, solution, witness
construction and independent checking. The stronger adaptive Taylor baseline
then checked every order and retained recurrence work. Its aggregate elapsed
time was 2.36 seconds versus 11.49 seconds for the earlier coarse-grid Taylor
scan, but these were initial, separately timed development runs with differing
success counts. This is baseline engineering, not a statistically established
speedup or learned-method result.

Operator squaring tightened 20/24 small-system residual bounds and could lower
Taylor order in two. A complete-cost test then used seven randomized paired
repeats on all 36 development cases, with the same checker in both arms. In the
eligible order-saving case, full squaring reduced degree 12 to 11 but increased
median complete runtime by about 11.3%. Squaring four selected anticommuting
groups reduced pair tests from 435 to 190; it remained about 2.3% slower in that
case. Both tested policies fail the efficiency gate.

A subsequent adaptive diagnostic reduced the block to 15 terms (105 pair
tests). Its approximately 0.6% median advantage occurred in only 10 of 15
paired repeats on the already studied case. This does not establish reliable
headroom or transfer. It is recorded as unresolved, not promoted to acquisition.

An additional bounded row-sum norm diagnostic examined 45 residuals using explicit small Hilbert-space matrices. It produced no new below-tolerance cases; its exponential construction cost is recorded and no larger-system scaling claim follows. Overlapping coefficient-split group certificates remain a mathematical candidate with no realized cost advantage.

The original grouping probe and first dense reference test contained defects;
they were explicitly invalidated and replaced. Their earlier apparent wins
are not counted. The mathematical checking and later timing results above use
the corrected implementations.

Evidence: [initial calculations](headroom_development.json),
[strong baseline](strong_baseline.json), [pruning](pruning_probe.json),
[full-square timing](quadratic_cost.json), [partial-square timing](block_cost.json),
[triple-block diagnostic](triple_block_diagnostic.json),
[protocols](../../research/v7/HEADROOM_PROTOCOL.md).

## Verification, costs and remaining obligations

The replay script exports exact polynomials and witnesses, reparses the saved
bytes and independently accepts 125 certificates. It also reproduces ten
strong-baseline budget refusals and the measured-data ambiguity witness.
Twenty-two V5/V6 source, result and data hashes match their prior receipts.
V4's calibration remains outside V5's dependency ledger.

Complete elapsed calculation times include failures and retries. Structural
operation counters and polynomial sizes are also retained, but are not a
complete bit-complexity or resident-memory model. Agent development time,
reference calculations and validation work are separate research costs; no
acquisition-inclusive net benefit is claimed. No new physical measurement
or hardware validation occurred.

[Exact certificate archive](certificate_archive.json) ·
[verification receipt](verification_receipt.json) · [tests](tests.log).

The goal remains active. Before a learner is built, a feasible construction
method must beat the strengthened baseline with sufficient complete-cost
margin across development problems. Only then may m1 be acquired and frozen,
followed by matched m2 acquisition with enabled, disabled, reconstructed,
restored and irrelevant controls; two untouched incremental improvements;
and a positive acquisition-inclusive net benefit. None of those learning gates
has been passed by the present results.
