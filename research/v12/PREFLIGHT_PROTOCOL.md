# V12 complete-cost preflight

Frozen before the first complete-cost paired run. The norm-only and checker-only
development diagnostics are already public. No held-out cases are used here.

Test the existing n=3 XXZ case, central Z, gamma=2, T=1/2, epsilon=1/1000,
512 live terms/cached columns and maximum degree 24. The baseline is V11
`guarded_taylor(..., guard='frobenius')`. The candidate is the same recurrence
with the V12 one-sweep cover hook using the union of ordinary weighted and
firstfit partitions and common denominator 2^20. A disabled-hook control
checks adapter overhead. The original V7 and V8 checkers remain unchanged.

Count model construction, every failed order/grouping/cover, accepted
construction, fresh complete checking and exact endpoint output. Record one
initial call per arm separately, then 21 three-arm blocks in shuffled order
using seed 1201. Preserve certificates, outputs, all timings and source hashes.
Output polynomials need not equal between different degrees. Disabled and
baseline polynomials, ordinary witnesses and exact endpoints must equal.

A useful preflight requires certified output plus at least a 5% reduction in
aggregate and median paired complete cost and wins in at least 18 of 21 blocks
against the strongest baseline. This is a development screening rule, not a
population-confidence theorem or a compounding claim. If it passes, prospectively
expand to the existing 36 development cases, counting failures and preserving
the same assumptions. If it fails, close this implementation; do not enlarge a
benchmark to rescue it or start acquisition. Timing diagnostics may identify a
separately justified mechanism, with a fresh declared preflight.

All algorithms and identities in this test are supplied conventional methods.
The reserved evaluation distribution remains unopened.
