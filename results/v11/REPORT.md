# V11: a cheaper conventional baseline, still no acquired method

The normalized-Frobenius rejection guard cuts total development calculation
cost by **36.8%**, including failures, while preserving every accepted
polynomial, norm witness, bound and endpoint. This is a supplied conventional
optimization. Compounding scientific capability remains unproved; neither
method acquisition nor reserved evaluation has run.

## Complete-cost comparison

The same36 inputs, physical assumptions, tolerance.001, degree24 limit and
original V7 checker were used in five repetitions. Each four-arm comparison
shuffled execution order. Construction includes failed order attempts and norm
proposals; fresh checking and exact endpoint evaluation are included.

| Method | Certified per36 cases | Total seconds over five repetitions | Successful paired wins over V8 | Median successful speedup |
| --- | ---: | ---: | ---: | ---: |
| Frozen V8 fraction-free Taylor | 26 | 10.940 | — | — |
| Max-coefficient rejection | 26 | 7.603 | 128/130 | 1.265× |
| Frobenius rejection | 26 | 6.912 | 130/130 | 1.346× |
| Max then Frobenius | 26 | 6.961 | 130/130 | 1.358× |

The simple Frobenius guard is the new conventional baseline. The two strongest
guards are close; these data do not establish an optimal guard schedule.
All ten failure cases remain failures. Across five repetitions each arm has
45 output-support refusals and five cached-column refusals; their cost stays
in the totals. No physical instance was made easier or removed.

## Mechanism and proof

For Hermitian R=sum r_P P on d-dimensional Hilbert space, Pauli orthogonality
and the singular-value bound give

`sum r_P² = trace(R²)/d <= ||R||infinity²`.

If the weighted Frobenius lower bound already exceeds the tolerance, no valid
upper norm certificate at that Taylor order can pass. The guard therefore
skips l1/firstfit/weighted proposals only at impossible orders. Equality
continues to the original proposals. This changes neither the first accepted
order nor the selected norm witness.

The implementation makes this decision directly from the common integer
derivative numerator using exact comparisons, before converting futile
derivative coefficients to Fractions. It still constructs all retained
polynomial coefficients and pays for all generator actions. On accepted runs,
the Frobenius guard skips910 orders across the campaign. Their construction
cost falls from4.530 to2.334 seconds; checking remains approximately2.45 seconds.

The lower bound is established mathematics, and applying it here is not
scientific novelty or autonomous m1 acquisition. Future comparisons must give
the strengthened baseline this operation. The [proof and implementation
audit](../../research/v11/lower_bound_gate.md) includes strict comparisons,
zero tails, derivative denominators and bounded integer temporaries.

## The matrix candidate lost

An exact d-by-d matrix action and Walsh-Hadamard Pauli conversion reproduced
the V8 polynomials, including nontrivial Y phases and local damping. Every
int64 kernel has a forward overflow bound and an exact Python-integer fallback.
No Liouville matrix or floating approximation was used.

On the prospectively fixed n3 XXZ, gamma2,T.5 case, it lost all seven complete
timing pairs. Mean cost was15.754ms versus13.839ms; the median ratio was1.133.
Construction and conversion overhead consumed the potential benefit. This
implementation is closed for the tested workload; the result is not a theorem
against matrix methods generally. Warmup costs are separately recorded.

## Verification and preserved evidence

- 720 guard-study calls, plus the seven matrix pairs and separately recorded
  initial calls.
- **106 saved certificates and endpoint outputs replayed** with the unchanged
  checker; exact cross-arm outputs match on all26 accepted physical cases.
- **190 tests passed**, including independent dense physical action, every
  two-qubit Pauli column, integer-overflow fallback, strict guard boundaries,
  rational comparison equivalence and refusal paths.
- Forty historical hashes, V7/V8/V9/V10 archives, and V8/V9/V10 source receipts
  remain unchanged. V4 calibration is still excluded from V5's dependency ledger.
- An initial agent profile measured the wrong fractional-cover path; it was
  rejected, preserved under `research/v11/invalidated/`, and replaced with a
  profile of the actual adaptive fraction-free baseline.

[Guard protocol](../../research/v11/REJECTION_GATE_PROTOCOL.md) ·
[Guard measurements](gate_cost.json) · [Matrix measurements](preflight.json) ·
[Replay receipt](verification_receipt.json) ·
[Gate archive](gate_certificate_archive.json).

## Next bounded method-construction test

Root analysis found that the old near-threshold cover does not need its full
69-group family or16 optimizer sweeps: the union of two ordinary partitions
and one sweep produces a checked14-group bound0.000998348 on the existing
degree11 residual. This refutes the initial negative motif inference, but its
norm-only timing is not end-to-end headroom.

The next test will use a common dyadic denominator for those split weights,
avoiding the old expensive rational reconciliation, then measure the full
calculation against the new Frobenius-guarded baseline. The existing V8 cover
checker remains fixed. See the [specific construction and stop
condition](../../research/v11/NEXT_COVER_PREFLIGHT.md). A narrower norm or cheaper
isolated subroutine still cannot substitute for method acquisition, causal
m1-to-m2 dependence, two untouched gains and acquisition-inclusive net benefit.
