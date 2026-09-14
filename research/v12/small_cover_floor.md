# V12 small-cover checker floor

This is a separate timing run for the actual `small_cover_taylor` implementation on the already-public case `n=3`, XXZ, `gamma=2`, `T=1/2`, `epsilon=1/1000`. The exact V11 baseline was produced with `guarded_taylor(..., max_order=12, guard='frobenius')`, and is degree 12. The new arm was produced by `small_cover_taylor(..., max_order=12, guard='frobenius')`; it certified at degree 11 with the V8 evolution-cover schema and exactly 14 groups (60 incidences). Both certificates passed their schema-selected independent checkers before timing.

Seven fresh paired checker plus exact endpoint operations gave these means in seconds:

| operation | V11 guarded degree 12 | actual small-cover degree 11 | saved |
|---|---:|---:|---:|
| checker | 0.0052689881 | 0.0051458156 | 0.0001231726 |
| exact endpoint | 0.0005404761 | 0.0004961544 | 0.0000443217 |
| checker + endpoint | 0.0058094643 | 0.0056420 | 0.0001674943 |

The small-cover arm saves about 0.123 ms in checking and 0.044 ms in endpoint evaluation, or 0.167 ms combined per paired run. These are checker-floor timings, not headroom claims. Small-cover construction time is recorded but excluded from this paired floor comparison.

The prior 61-group diagnostic is retained only as `oldwitness` historical context. It was not used for ranking, cost comparison, or witness selection. Full per-repetition receipts, actual group and incidence counts, certification status, construction breakdown, and source hashes are in [small_cover_floor.json](/Users/aidenlippert/Documents/Spectra/results/v12/small_cover_floor.json).
