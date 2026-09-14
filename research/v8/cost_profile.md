# V8 cost profile of complete adaptive Taylor certificates

This is a read-only development profile of `experiments/v7_adaptive_taylor.py`
and `experiments/v7_certificate.py`, using the cached Python executable named in
the README and NumPy-only project dependencies. Each case included proposal
generation, witness construction inside the adaptive proposer, and a fresh
`Generator` plus `check_certificate` replay. `cProfile` and `tracemalloc` were
used; no existing output writer was rerun or modified.

## Four successful cases

Wall times include profiler overhead, so the phase ratios are the useful
comparison. Counters are exact source-level work counters from `Generator` and
the independent checker.

| n, gamma, T | order / polynomial entries | generation s | checker s | peak bytes | proposal columns / cache hits | coefficient adds | peak live terms | checker pair checks / squares |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 3, 0, 1/5 | 5 / 53 | .0284 | .0222 | 111,342 | 53 / 23 | 210 | 18 | 0 / 36 |
| 3, 2, 1/2 | 12 / 303 | .1128 | .0656 | 172,274 | 303 / 273 | 1,449 | 30 | 0 / 60 |
| 4, 0, 1/5 | 5 / 117 | .0888 | .0741 | 230,725 | 117 / 30 | 631 | 65 | 0 / 130 |
| 4, 2, 1/2 | 13 / 1,175 | .6499 | .3214 | 635,556 | 1,175 / 1,049 | 7,219 | 126 | 222 / 158 |

The complete profiled elapsed times were .0506, .1783, .1629 and .9714 s,
respectively. Generation was 56%, 63%, 55% and 67% of those totals; checking
was the remaining 44%, 37%, 45% and 33%. Thus norm tightening cannot be called
the dominant opportunity for the whole pipeline from these cases: in the
largest case `apply` plus exact arithmetic dominates proposal time, and fresh
residual recomputation dominates checking time.

The cumulative profile identifies the concrete kernels. In the largest case,
`adaptive_taylor` took .650 s; its 28 `Generator.apply` calls took .501 s.
`check_certificate` took .321 s, with `residual_records` at .315 s and
`norm_witness` at .313 s. Fraction construction/add/multiply accounted for
substantial nested time (.265/.155/.147 s cumulative for the most visible
operations), while `_anti` took .170 s and Pauli commutators .125 s. The
smaller gamma=0 cases likewise put `apply`/`commutator_i` and exact arithmetic
ahead of grouping. Norm grouping becomes material at n=4, gamma=2, but it is
one part of checker replay rather than the universal bottleneck.

## n6 refusal and representation lifetime

At the v7 cap (`max_terms=512`), the observed adaptive outcomes were:

| n=6 case | first refusal | counters at refusal |
|---|---|---|
| gamma=0, T=1/5 | succeeds | — |
| gamma=0, T=1/2 | `cached column budget` | cache 512, peak live terms 456, 739 requests, 226 cache hits |
| gamma=2, T=1/5 | `output term budget` | cache 340, peak live terms 425, 621 requests, 281 hits |
| gamma=2, T=1/2 | `output term budget` | same first refusal counters as gamma=2,T=1/5 |

The gamma=0, T=1/2 refusal is therefore partly a representation/lifetime
artifact: increasing the cap to 1024 still refused at cache 1024 (peak live
terms 707), but cap 2048 completed at order 11 with 4,560 polynomial entries,
2,009 cached columns and peak live terms 1,053. The cache is retained across
all adaptive norm probes and orders, so its lifetime is longer than the live
support needed by any one operator. This does not prove that eviction is safe:
it identifies a specific reusable-column lifetime seam to test.

For gamma=2, T=1/5, caps 1024 and 2048 show the other seam: cap 512 refuses
with 340 cached columns and 425 live output terms; cap 1024 refuses with 678
cached columns and 741 live output terms; cap 2048 succeeds at order 7 with
1,527 polynomial entries, 741 cached columns and peak live terms 1,163.
For gamma=2, T=1/2, cap 2048 succeeds at order 16 (18,355 polynomial entries,
2,046 cached columns, peak live terms 2,046). These cap experiments were
bounded diagnostic reruns and did not write project outputs.

## Specific mathematical/algorithmic headroom

1. **Separate cache lifetime from live-support budget.** The current cache is a
   monotone set of columns, while `apply` also enforces the instantaneous output
   support cap. A sound generation schedule could retain columns by future-use
   count or recompute cold columns, with a proof that the resulting exact
   certificate is unchanged. The gamma=0,T=1/2 transition (512 -> 2048) makes
   this a concrete target; it is not a generic caching suggestion.

2. **Avoid recomputing residual records for an independent checker.** The
   checker intentionally recomputes every residual, which costs .315 s in the
   largest case and repeats the same 28 generator applications. A compact
   algebraic transcript (operator hashes plus exact recurrence checkpoints)
   could let the checker verify the proposer’s recurrence while preserving an
   independent recomputation of selected checkpoints. This requires a new
   soundness argument; it is a certificate-design opportunity, not merely a
   faster implementation.

3. **Target commutator closure and rational growth.** The dominant proposal
   work scales with coefficient multiply/adds (7,219 in the largest case),
   exact Fraction operations, and Pauli commutator construction. A mathematically
   sparse closure representation that tracks only nonzero reachable Pauli labels
   with exact integer numerators and shared denominators could reduce repeated
   normalization while preserving exactness. Standard integer scaling or
   caching alone is not treated as learned novelty.

Norm tightening remains a valid local experiment: `norm_witness` is visible in
the n4,gamma=2 checker profile and `_anti` costs .170 s cumulatively. But the
complete profiles do not support making it the dominant opportunity without
also accounting for residual recomputation, commutators, and rational
arithmetic.

## Root audit correction

The listed `norm_witness` .313 s attribution to checker replay is invalid: the checker never calls that function. Treat it as an unseparated profile statistic, not a checker subphase. The raw four-case profile still identifies exact generator and Fraction work, but detailed additive attribution is not established. Sampling only selected recurrence checkpoints would require a probabilistic soundness theorem and cannot silently replace deterministic residual checking. No such checker change was made. Cache eviction was subsequently tested at unchanged live support and recovered zero cases; increasing support to 2048 is not the same budget.
