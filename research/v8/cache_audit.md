# V8 bounded-cache audit

This audit tests whether evicting cached generator columns can recover any of
the 36 fixed v7 development cases while preserving the v7 requirements:
`max_terms=512` for instantaneous operator support, polynomial degree at most
24, tolerance `1/1000`, and the original six `(gamma,T)` settings. The
Hamiltonian families and widths are the same `v7_headroom.model` cases. No
core files were edited.

The temporary subclass used a 512-entry LRU column cache. On a miss at a full
cache it evicted one column and recomputed it exactly with the existing
`commutator_i`, `Fraction`, and output-support checks. Hits were refreshed in
LRU order. The checker also used this subclass. Thus cache capacity stayed
512; recomputation was charged in the existing source counters, and termcap
was never increased or treated as equivalent to cache capacity.

Result: **zero of 36 cases were recovered.** All 24 width-3/4 cases and the two already successful width-6 cases
retained their baseline success status; the ten refused cases remained refused. The
cache seam is therefore not sufficient at the fixed live-support budget.

The one direct cache-pressure target, width 6 XXZ, `gamma=0`, `T=1/2`, did
evict 88 columns, but subsequently refused on `output term budget` with peak
live support 456 during the adaptive run. The original retained-cache run
refused earlier on `cached column budget` (cache 512, peak live support 456).
Eviction changes the refusal reason and permits more requests; it does not
reach an accepted degree-24-or-less certificate. The later live-support
growth is the binding seam. Other width-6 refusals already hit output support
before cache capacity, so eviction cannot help them.

The audit's temporary LRU policy is not a proof that no more sophisticated
future-use policy can help. It does establish that bounded eviction, as a
straightforward remedy for the observed retained-cache lifetime, does not
recover a fixed case under the same live support cap. Any claim of recovery
would need a new schedule and a separate exact checker replay, with its
recomputation cost reported.

Profile claim check: `check_certificate` calls `residual_records` and validates
the supplied witness groups; it does **not** call `norm_witness`. Therefore a
statement that checker replay contains `norm_witness` at `.313 s` as a nested
measurement is not supported by this call graph. `norm_witness` is called by
the adaptive proposer (and by `derive_certificate`), while
`residual_records` is called by both proposer-side certificate construction
and checker replay. The `.313 s` figure should be described as a cumulative
profile statistic unless the omitted profiling harness demonstrates a
separate checker invocation that calls it. It cannot be a direct nested
checker call in the code audited here.

Exact rows and counters are in [`results/v8/cache_audit.json`](../../results/v8/cache_audit.json).
