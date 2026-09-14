# V10 headroom audit: exact-solvable-reference Duhamel method

## Scope and evidence

This is a bounded read-only preflight. It uses the already fixed 36 physical
development cases, tolerance `1/1000`, horizons `T=.2/.5`, live support/cache
cap 512, and degree cap 24. It does not open or benchmark reserved n5/n7
workloads. The comparison target is the fraction-free adaptive Taylor/BFS
family already measured in V9. The inputs are [V7's complete-cost baseline](../../results/v7/REPORT.md),
[V8's exact norm and cost replays](../../results/v8/REPORT.md),
[V9's trajectory/collocation replay](../../results/v9/REPORT.md),
[the reference-semigroup derivation](../v9/reference_semigroup.md), and the
[fixed V7 workload driver](../../experiments/v7_headroom.py).

The strongest relevant timing evidence is V9: fraction-free adaptive Taylor
certified 26/36 at 2.176 seconds; BFS certified 27/36 at 4.240 seconds;
trajectory enrichment was 27/36 at 7.769 seconds. The 26 baseline successes
and ten refusals therefore define a difficult, already measured acceptance
frontier. V7/V8 also show that smaller degree, norm tightening, overlap covers,
cache policies, and incremental reuse did not produce complete-cost headroom.

## What the factorial bound buys

For `G=G0+B`, the reference construction has residual `-B D_m` and certified
tail

```
             ||O|| (K T)^(m+1)/(m+1)!,   K >= ||B||_(infinity -> infinity).
```

This is useful at `T=.2` when `K` is genuinely small. The supplied reference
note gives the actual n6 XXZ triangle bound `K <= 8.8`, hence `K T <= 1.76`.
The scalar tail at orders 4, 6, 8, and 10 is respectively approximately
`1.41e-1`, `1.04e-2`, `4.47e-4`, and `1.26e-5`. Thus order 8 clears `1e-3`
under this bound (before any initial, conversion, or certificate margins).
That is the only plausible headroom region: short horizon and a reference
whose interaction norm is substantially smaller than the full generator.

At `T=.5`, the same bound has `K T <= 4.4`; its order-8 tail is about 1.703,
and order 12 is about 0.03721. The factorial estimate cannot certify the
requested tolerance there at degree 12. Its first passing order is 15
(`4.4^16/16!` is below `999/1000000`), below the degree24 budget. The earlier
agent estimate of order19 and proximity to the cap was incorrect.
Mixed cases retain additional onsite/YY interaction in `B`, so the XXZ number
is not a favorable-family-wide estimate. At `T=.2`, a modest `K T=2.4`
already needs order 10 just to get below `4e-4`; at `K T=3.0`, order 12 is
about `2.6e-4`. These are scalar bounds, not measured errors, and cannot be
used to claim success.

## Why complete-cost headroom is unlikely

The exact reference removes diagonal evolution cost but does not remove the
dominant combinatorics. Each Duhamel level applies `B` to every retained
operator mode and then convolves it with every `G0` frequency. The number of
terms is driven by Pauli support growth and frequency collisions; collisions
also create polynomial factors. The fixed support cap 512 turns this into
refusal work rather than free truncation. A sound implementation must charge
mode generation, repeated-frequency handling, exact coefficient growth,
outward exponential/trigonometric enclosures, failed orders, conversion to
the existing checker, and independent replay.

The reference therefore trades Taylor recurrence work for a larger symbolic
object. V9 already measured that a trajectory basis with 213 labels can be
slower than a 181-label construction (BFS .147 s versus enrichment .281 s on
the shared n6 short-horizon case), and that a lower-degree collocation object
still lost to Taylor because solve/conversion/checking dominated. Duhamel
mode count is at least as exposed to this effect, while its useful factorial
tail gets more expensive at `T=.5`. Any complete-cost claim must use paired
contemporary timings and include failures. Beating a historical wall-time
total from a different run is not the gate. The preflight establishes no
cost advantage; this is a prediction, not an impossibility theorem.

## One concrete alternative if the preflight fails

Use a **two-layer commutator-free exponential stepper** instead of expanding
all Duhamel modes. Keep the exactly diagonal `G0`, and partition the supplied
interaction into disjoint even and odd edge layers, grouping the commuting XX
and ZZ terms on each edge (and the analogous terms in mixed cases). On each
step apply the exact local/tensor-product exponentials for

```
S_h = exp(h G0/2) exp(h B_even/2) exp(h B_odd)
      exp(h B_even/2) exp(h G0/2).
```

The construction retains a fixed sparse operator map rather than a growing
exponential-polynomial mode list. Its certificate can use the existing
variation-of-constants checker with the exact defect
`S_h' - G S_h` represented by commutator terms, or a separately proved
outward commutator envelope. This only merits pursuit if that envelope is
computed from the actual Pauli commutators and all step/refusal work is
charged; generic scalar splitting bounds are insufficient. It is a concrete
representation change because each layer has O(n) local factors and does not
enumerate Duhamel frequency histories.

The smallest decisive test is one existing case: n3 XXZ, `gamma=2`, `T=.5`,
with the same tolerance, support/cache and independent checker. Run the
stepper over step counts `h=.1,.05,.025`, include every failed step and
commutator-bound construction in elapsed cost, and require a certified result
whose complete time is below the matched fraction-free Taylor call. If none
of the three steps both certifies and wins, stop the alternative; do not
broaden it to reserved sizes. This test is deliberately on the fixed physical
workload and does not assume favorable noise or use held-out data.

## Decision

Treat exact-reference Duhamel as a short-horizon diagnostic only. Before any
learner acquisition, require the full 36-case replay with construction,
mode growth, exact enclosures, failures, conversion, and independent checks
included. On present evidence its factorial advantage is likely consumed by
mode/coefficient explosion at the long horizon and by symbolic overhead even
where the short-horizon scalar tail is small. If the smallest split-step test
fails the complete-cost gate, this reference line should be closed rather
than rescued by a looser norm, more agents, or a renamed known method.
