# V11 fraction-free Taylor baseline

Authoritative workload: existing `v7_headroom` model, `n=3`, `family=xxz`,
`gamma=2`, `T=1/2`, `epsilon=1/1000`. The measured proposer is exactly
`experiments.v8_integer_taylor.fraction_free_taylor`, at order 12. Each run
then constructs a fresh original V7 `Generator`, calls `check_certificate`
with the returned piece and witness, and evaluates `p.coefficients` exactly at
`T`. No fractional-cover proposer, held-out data, tuning, or baseline edits
were used.

There was one warmup and seven paired repeats in one process. Mean timings:

| phase | mean seconds | median seconds |
|---|---:|---:|
| fraction-free construction | 0.007852 | 0.007820 |
| fresh V7 check | 0.005516 | 0.005487 |
| exact endpoint evaluation | 0.000530 | 0.000529 |
| all components | 0.013898 | 0.013859 |

All seven checks were `certified`, with bound
`9235404595038841/20270250000000000000`. The exact endpoint coefficient for
the initial observable and all operation counters, source hashes, raw rows,
and one cProfile trace are recorded in
[`baseline_profile.json`](../../results/v11/baseline_profile.json).

The construction phase is 56.5% of the all-component mean; fresh checking is
39.7% and endpoint evaluation is 3.8%. If proposer construction were made
free while the checker and endpoint remained unchanged, the measured ceiling
would be `0.013898 / 0.006046 = 2.30x`, or at most 56.5% time reduction for
this workload. This is a ceiling, not an achieved speedup. Counters in the
JSON are cumulative operation counts; phase seconds are wall-clock self-time.
The cProfile `cumtime` column is inclusive call time and `tottime` is function
self-time.

The earlier fractional-cover measurement was invalid for this question and is
preserved under `research/v11/invalidated/` for auditability.
