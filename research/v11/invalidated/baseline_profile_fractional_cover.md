# V11 baseline profile: proposer-only ceiling

This is a bounded contemporaneous baseline for the existing `v7_headroom`
model (`n=3`, `family=xxz`, `gamma=2`, `T=1/2`, `epsilon=1/1000`). It uses
the complete order-8 V8 fractional-cover calculation on the V7 full Taylor
coefficients, followed by a freshly constructed V7 witness, an independent
fresh V7 check, and exact rational endpoint evaluation. No held-out data,
tuning, or baseline code changes were used.

The measurement used the bundled Python runtime, one warmup, then seven
paired repeats in one process. Each repeat ran the same four phases in order.
`*_seconds` are phase wall-clock self-times; `total_seconds` is their sum.
The operation counters are cumulative counts from the relevant generator or
checker, not self-time measurements. A cProfile trace was collected for one
additional complete run; its `cumtime` is inclusive call time and `tottime`
is function self-time. Raw rows, counters, and the profile excerpt are in
[`baseline_profile.json`](../../results/v11/baseline_profile.json).

The seven-repeat means were:

| phase | mean seconds | median seconds |
|---|---:|---:|
| V8 construction + proposer | 0.09010 | 0.08366 |
| V7 witness construction | 0.00451 | 0.00390 |
| fresh V7 check | 0.00532 | 0.00385 |
| exact endpoint evaluate | 0.00040 | 0.00033 |
| all measured phases | 0.10033 | 0.09216 |

The endpoint coefficient for the initial observable was the exact rational
`1001445889/3937500000`. The V7 check returned `over_tolerance` in all seven
repeats, with reported bound
`1902120580047273743/46080000000000000000`; this is not a certified green
check and is not treated as learning success.

The construction phase is about 89.8% of the measured mean. If a hypothetical
proposer acceleration removed the entire measured construction phase while
leaving witness, checking, and endpoint costs unchanged, the arithmetic upper
ceiling is `0.10033 / 0.01023 = 9.81x` (about 89.8% maximum time reduction).
This is a ceiling for accelerating the proposer portion of this workload, not
an achieved speedup and not evidence that the remaining phases can be removed.
The corresponding non-construction floor is about 10.2 ms per repeat.

The cProfile run attributes the proposer’s inclusive time primarily to
`propose_fractional_cover` (0.218 s), `_group_family` (0.046 s),
`rationalized` (0.065 s), and the proposer’s internal exact
`check_fractional_cover` calls (0.040 s). This makes the proposer the clear
measured bottleneck, while also showing that its internal exact verification is
already charged to construction. The profile is a diagnostic trace; the
seven-repeat paired wall-clock statistics are the timing claim.
