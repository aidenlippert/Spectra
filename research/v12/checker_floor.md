# V12 checker floor

This is a certified, already-public case: `n=3`, XXZ, `gamma=2`, `T=1/2`, and `epsilon=1/1000`. The degree-12 baseline was produced by the exact call

```python
guarded_taylor(Generator(h, gamma, 3, 512), initial, T, TOL,
               max_order=12, guard='frobenius')
```

where `h, initial = model(3, 'xxz')`. It returned degree 12, with claimed integrated bound `9235404595038841/20270250000000000000`. Its coefficient support counts by degree were `[1,3,9,20,30,30,30,30,30,30,30,30,30]`. The measured construction total was recorded separately from checker and endpoint timing; the detailed integer and guard work is in `results/v12/checker_floor.json`.

The accepted V8 overlap arm uses `p11 = Piece(T, p12.coefficients[:12])`. The final residual was independently derived through `residual_records` (equivalently `-G c11`). It has 30 terms. The existing passing row 10 in `results/v8/fractional_probe.json` has exactly the negated residual, so its cover coefficients were negated before use. The exact evolution-cover checker certified the resulting witness with integrated bound `19761823751632003/20480000000000000000`.

Seven fresh paired repetitions measured checker plus exact endpoint evaluation. Means (seconds) were:

| operation | degree-12 guarded baseline | degree-11 V8 overlap | baseline minus overlap |
|---|---:|---:|---:|
| checker | 0.0053078631 | 0.0075265891 | -0.0022187260 |
| exact endpoint | 0.0005524643 | 0.0005089940 | +0.0000434703 |
| checker + endpoint | 0.0058603274 | 0.0080355831 | -0.0021752557 |

On this fresh local measurement, replacing the degree-12 guarded checker with the degree-11 overlap checker saves no time: the overlap checker costs about 2.22 ms more on average. Endpoint evaluation saves about 0.043 ms, leaving a combined increase of about 2.18 ms. This is a checker-floor result only; it is not headroom. The overlap proposer’s construction cost is excluded from this floor timing, as requested, and no held-out calculation is used.

The complete machine-readable receipt, including all seven rows, checker costs, construction breakdown, exact calls, and source hashes, is [checker_floor.json](/Users/aidenlippert/Documents/Spectra/results/v12/checker_floor.json).
