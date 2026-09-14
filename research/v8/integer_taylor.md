# Fraction free Taylor baseline

Implemented `experiments/v8_integer_taylor.py` as a conventional exact
baseline. It clears the common denominator `d` of the Hamiltonian generator,
clears the initial denominator `q0`, and propagates integer sparse numerators
with `b[k+1] = (dG)b[k]`. Exported coefficients are exactly
`c[k] = b[k] / (q0*d**k*k!)` as `Fraction` values. The proposer uses the same
exact residual norm witness shape as v7; the unchanged independent
`check_certificate` is run with a fresh `Generator` in the test.

The implementation keeps the existing order and support/cache budgets at 24
and 512. It reports denominator and integer conversion counts in its
construction cost, while generator construction, witness grouping, and
checking remain explicit separate phases for callers to time. This is a
conventional denominator-clearing baseline only; it is not claimed as an
acquired research operation or m1 novelty.

The recurrence now caches cleared integer columns lazily (bounded by the same
512 cache cap), propagates `b` with Python integers, checks instantaneous output
support, and converts only the resulting derivative to exact fractions for the
existing norm witness. It does not call `gen.apply(c)` in the integer hot loop.
Input operators are passed through the existing `clean` validation and cleared
denominator bit length is capped at 8192 bits.

Validation command (cached README Python):

```text
/Users/aidenlippert/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 -m unittest tests/test_v8_integer_taylor.py -v
```

Result: `Ran 4 tests ... OK`. The tests compare exact `Piece` coefficients and
complete claimed bounds against `adaptive_taylor` for development cases
`n=3, gamma=0, T=1/5` and `n=4, gamma=2, T=1/2`, then independently certifies
both generated witnesses with a fresh v7 checker generator, plus stationary
rational initial data, cache/live refusal behavior, and invalid labels.
