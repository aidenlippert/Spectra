# Better selection, with a first collective ceiling for the completed family

**Trace-normalized ranking reduced the certified H6 interval to 3.415421 mHa,
versus 3.942479 mHa for the preceding coefficient-ranked run.** Both started
from the same 64-combination certificate under a 360-second discovery budget.
The final counts were 216 and 214, respectively. At a common ceiling of 214
combinations, the new run's 208-combination certificate gives **3.435563 mHa**.
The improvement therefore persists without allowing the new run more
combinations.

The conditional diagnostic also produced the **first accepted exact dual
witness covering the entire 492-generator completed-spin family**, including
every unselected direction. With the frozen reference upper, its ceiling
implies that this family cannot produce an interval narrower than
**1.513465 mHa**. This does **not** rule out the 1.6 mHa target, which the
discovery run still missed.

The original compression question remains open. The ten spatial Hamiltonian
patterns and their **0.108331264 mHa collective remainder bound** stayed fixed.
The growing objects are correlation combinations and their proof coefficients.
216 combinations are 43.90% of this 492-generator frame; that is not yet
evidence that a small representation captures all difficult correlations.

![Comparison of the two ranking rules at common resource ceilings](/Users/aidenlippert/Documents/Spectra/results/trace_pricing_20260913/comparison.png)

## What the comparison establishes

The old ranking used the Euclidean size of a direction's frame coefficients.
The new rule divides its proposed moment violation by its mean squared
operator size under the uniform fixed-N trace. Both rules then test the
actual energy gain through the same solver, rational export and exact checker.
The negative eligibility threshold is interpreted in the new normalized
units; this pass tests that whole proposal rule, not the metric in isolation.

The controller, two candidate bundle sizes, gain-per-cost selection,
60-second reserve, 40-second solver limit and rounding denominator were
unchanged. The controller's only code difference is the imported pricing
model. All nineteen new candidates passed fresh exact replay. All nine
completed two-bundle comparisons selected the larger bundle, so this pass
does not demonstrate an advantage for adaptive bundle size itself.

| Combination ceiling | Old best interval, mHa | New best interval, mHa | New count used |
|---|---:|---:|---:|
| 64 | 7.414219 | 7.414219 | 64 |
| 96 | 6.193385 | 6.006426 | 96 |
| 128 | 5.240350 | 4.792560 | 128 |
| 160 | 4.586025 | 4.068281 | 160 |
| 192 | 4.186035 | 3.649392 | 192 |
| 214 | 3.942479 | 3.435563 | 208 |
| 256 | 3.942479 | 3.415421 | 216 |

| Incremental time ceiling | Old best interval, mHa | New best interval, mHa |
|---|---:|---:|
| 60 s | 5.804984 | 5.630248 |
| 120 s | 4.672596 | 4.222983 |
| 180 s | 4.312888 | 3.807123 |
| 240 s | 4.111738 | 3.536212 |
| 300 s | 4.087226 | 3.435563 |
| 360 s | 3.942479 | 3.415421 |

These are the best certificates actually completed within each ceiling;
unrun results are not interpolated. The common source appears at zero
incremental time, and its prior discovery cost is charged separately below.
This is one historical comparison, with the old run frozen rather than
repeated. It establishes a measured improvement here, not statistical
reproducibility or general superiority. The full tables preserve every
candidate, including those not selected.

The new final interval is 13.37% narrower than the old final interval and
53.93% narrower than the common source. The first sampled new result to beat
the old final interval used **168 combinations**, giving 3.939360 mHa at
140.512 seconds. The last step, from 208 to 216, gained only 0.020142 mHa.

The best lower is -6.336474047298517... Ha, paired with the same rational
upper -6.333058626233001... Ha. The numerical solver returned
`optimal_inaccurate`; the reported lower includes its exact coefficient-L1
residual penalty of **0.009140340 mHa**. Acceptance rests on the reconstructed
rational proof, not the numerical status.

## Cost and numerical stability

| Quantity | Common 64 seed | Old 214 result | New 216 result |
|---|---:|---:|---:|
| Gram entries, including quadratic baseline | 6,236 | 11,450 | 11,556 |
| Numerical map nonzeros | 251,268 | 2,537,322 | 2,579,604 |
| Compact lower certificate | 94,565 B | 237,366 B | 252,023 B |
| Expanded SOS certificate | 645,205 B | 1,453,772 B | 1,477,698 B |
| Nonzero direction coefficients | 3,936 | 13,224 | 13,284 |
| Nonzero compact factor coefficients | 4,736 | 9,402 | 9,606 |

The new selected blocks all have dimension 27. Its last trial cost 44.991
seconds: 1.273 construction, 36.482 solve, 0.034 export and 7.202 exact
acceptance. The full new campaign cost **319.227 seconds**, versus 332.067
for the old continuation. It stopped on the declared time reserve; the
256-combination and twelve-round caps were not reached. Counting the common
222.511 seconds of causal source discovery gives **541.738 seconds** for
the new branch and 554.578 for the old branch. Earlier failed controls and
validation costs remain in their original ledgers.

Trace-metric setup added **0.011262 seconds**, and all ten pricing calls
together took 0.103172 seconds. The raw coefficient trace metrics are poorly
conditioned, with a maximum condition estimate of 6.22e10. Direct generalized
eigensolves failed the coordinate-rescaling test; diagonal equilibration
alone also failed. Both failures and the intermediate source versions are
preserved. The final implementation constructs the pricing Grams in the
orthonormal operator-coefficient basis from a full thin QR decomposition.
Its maximum metric condition estimate is 116.804, and the unchanged
rescaling tolerance passes. No ridge, rank truncation or change to the
energy solver's conditioning was used.

Fresh energy replay and the control comparison took **145.250 seconds**.
The separately bounded full-family dual diagnostic cost **203.310 seconds**.
They ran concurrently after discovery; their wall times are not additive
calendar time. Neither competed with the timed discovery run. BLAS thread
environment variables were one; Clarabel retained its automatic setting.

## A certificate covering all unselected combinations

The final numerical dual was projected onto the exact affine constraints.
Uniform-trace mixing then supplied a proposed repair. The estimated minimum
mixture was 0.132500%; the first scheduled proposal, **0.135150%**, passed
all exact gates, including the eight complete anticommutator PSD blocks of
dimensions 30, 30, 93, 93, 93, 93, 30 and 30. Its witness is 86,963 bytes.
No full-family primal SDP was run.

The accepted ceiling on every molecular lower in this proof family is

    C = -6.334572090891981... Ha.

It applies to the full completed frame and all its selected subspaces,
with the quadratic baseline, body-one number ideal, coefficient-L1 residual
and fixed lower tail shift used here. It is not a ceiling on arbitrary
many-body methods or on every certificate built from the same ten spatial
patterns.

Since every such lower L obeys L <= C, the frozen upper U gives

    U - L >= U - C = 1.513464658980... mHa.

The current lower remains 1.901956 mHa below that ceiling. A successful
1.6 mHa interval would require a lower in the remaining **0.086535 mHa**
window immediately below this ceiling. The witness neither proves that
this window is attainable nor rules it out.

Independent standard-library replay of this witness and the prior precise
physical proof passed in **181.580 seconds**. The latter proves
E0 >= -6.333113620625836... Ha for the identical rational Hamiltonian and
sector. Therefore every lower L in this family has **E0 - L >= 1.458470 mHa**.
This physical-error floor differs from the fixed-upper interval floor; it
also does not exclude 1.6 mHa. The physical reference itself encloses E0 in
0.054994 mHa, which bounds the remaining uncertainty in the frozen upper.

## Next decisive experiment

Freeze this comparison and **extract a reusable generation rule from the
successful coupled operators**, as requested in the user's subsequent
steering. Diagnose lower-family, upper-reference, representation and exact
export contributions separately. Test the candidate structure on untouched
geometries and a larger fixture under bounded budgets, with the rule frozen
before inspecting their results. Reaching 1.6 mHa on this H6 instance is not
a prerequisite for testing transfer.

The present ceiling leaves family adequacy unresolved. A stronger repair
could settle it, but it will not displace the mechanism experiment with
another sequence of H6 enrichment runs. The current result does not justify
abandoning the family or declaring that more selected combinations will
succeed. No general compression theorem or scaling result follows from
this frozen finite-basis example.

## Verification and reproduction

Six focused tests pass: exact trace-metric agreement, coordinate-rescaling
invariance, rounded direction validity and independence, refusal of invalid
inputs, and two tests of actual count/time ceilings and incumbent retention.
The first numerical failures are preserved alongside the final passing logs.

All nineteen new lower certificates, the common 64-combination source, and
the frozen control's best lower passed fresh exact replay against the same
rational upper. Other historical checkpoints are tied to their previous
exact replays by matching hashes. The accepting comparison loaded no NumPy,
SciPy, CVXPY or PySCF. All **385 files** in the prior manifest are unchanged.
No many-body sector was enumerated for the new lower or dual proofs. The
reference upper still uses its previously discovered 200-amplitude witness;
that discovery has not been removed by this method.

The [protocol](/Users/aidenlippert/Documents/Spectra/research/trace_pricing_20260913/PROTOCOL.md)
and [proof](/Users/aidenlippert/Documents/Spectra/research/trace_pricing_20260913/PROOF.md)
specify the changes and limits. The
[fresh comparison](/Users/aidenlippert/Documents/Spectra/results/trace_pricing_20260913/fresh_replay.json),
[cost ledger](/Users/aidenlippert/Documents/Spectra/results/trace_pricing_20260913/cost_ledger.json),
[full-family witness](/Users/aidenlippert/Documents/Spectra/results/trace_pricing_20260913/full_dual/witness.json)
and [fresh physical-error audit](/Users/aidenlippert/Documents/Spectra/results/trace_pricing_20260913/audit.json)
hold the exact results and resource records.

From the workspace root, with the configured Python environment:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest \
  research.trace_pricing_20260913.test_pricing \
  research.trace_pricing_20260913.test_comparison -v
python -S -m research.trace_pricing_20260913.replay \
  --results results/trace_pricing_20260913 --out /tmp/trace-pricing-replay.json
python -S -m research.spin_enrichment_20260913.audit \
  --results results/trace_pricing_20260913 --out /tmp/trace-pricing-audit.json
```

`research.trace_pricing_20260913.campaign --out <new-directory>` repeats
the six-minute trace continuation. Existing discovery directories are
refused. The SVG beside the plotted PNG is available for export.
