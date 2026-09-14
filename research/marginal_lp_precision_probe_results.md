# LP-to-exact precision probe

This probe reruns the asymmetric polishing LP in an isolated script and measures where the HiGHS objective is lost during exact export. It does not alter the production polisher or reuse its output directory.

For `1_1000_penalty/polished/certificate.json` (253 candidate squares), HiGHS returned objective 3.2809863869603477 with maximum floating equation residual 1.0701e-7. The exported exact certificate replayed at 3.2809749582861105, a gap of 1.1428674e-5. No polynomial terms were dropped by the 1e-14 coefficient cutoff (4015 retained, 0 dropped), so that export cutoff is not the cause. Exact replay reports residual L1 1.1428358e-5, essentially the entire objective gap. The run used 15 orbit and 236 direct squares.

The same diagnostic on `1_1000_atom_exact/certificate.json` (1184 candidates) gives HiGHS objective 3.2805934702722936, max equation residual 9.5296e-8, exact lower 3.2805832671114383, and objective gap 1.0203161e-5. Again, zero polynomial terms were dropped (18,365 retained); exact residual L1 is 1.0203158e-5. The weighted equation L1 values are 1.1417134e-5 and 1.0202626e-5 respectively.

The assembled LP matrices contain extremely small nonzero entries: minimum absolute values are 3.36e-26 (penalty) and 7.45e-37 (atom), with 1% quantiles 9.60e-22 and 3.36e-18. A same native HiGHS 1.15.1 backend experiment at fixed 253 rays isolated the cause. With `small_matrix_value=1e-9`, the replayed lower was 3.28097496857942 and width 2.9809304e-5; with only that option changed to `1e-12`, the lower became 3.2809853575378245 and width 1.9420346e-5. The numerical-to-export gap fell from about 1.14e-5 to 9.01e-10, while replay residual L1 fell from 1.1418e-5 to 1.9532e-7. This supports internal small-entry dropping as the causal source of the earlier loss, rather than polynomial coefficient rounding. The option and its default are documented by HiGHS: [option definitions](https://ergo-code.github.io/HiGHS/dev/options/definitions/).

The earlier SciPy wrapper limitation remains accurate for SciPy 1.14.1's public API: it does not expose a documented `small_matrix_value` parameter. The native `highspy` interface does expose it and makes the controlled attribution possible.

The isolated outputs are in the per-source directories `results/marginal_lp_precision_probe/1_1000_penalty_polished/` and `results/marginal_lp_precision_probe/1_1000_atom_exact_certificate/`. The root-level probe files are historical overwritten diagnostics, not separate evidence for both sources. The controlled native comparison is in `results/marginal_native_precision/1e-09/` and `1e-12/`. Reproduce with:

```text
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_lp_precision_probe results/marginal_asymmetric_adapt/1_1000_penalty/polished/certificate.json
```

These measurements diagnose the gap; they do not certify the HiGHS numerical optimum. A certified bound must use the exact replayed certificate (or a separately verified rational repair).

The production polishing entry point now defaults to native HiGHS with `small_matrix_value=1e-12`. Rerunning all 1,196 raw adapted directions produced 13 orbit squares plus 237 direct squares, lower 3.280985482140053, upper 3.2810047778834033, width 1.929574335056606e-5, and residual L1 2.1192199212155313e-9. Its accepted artifacts are in `results/marginal_native_precision/polished/`; a separate `python -S` replay and exact comparison against the epsilon=1/1000 model passed.

A controlled weight-cap sweep on the same raw directions tested caps 100, 1,000, and 10,000. Exactly one ray reached the cap in each run. The numerical objective changed only about 1.7e-8 across the sweep, while exact widths worsened from 1.9295743e-5 to 1.9305307e-5 and 1.9402087e-5. Increasing the cap does not close the remaining gap on these fixed directions.

The saturated raw ray is nearly a fixed-number null polynomial. After correctly summing all CAR action paths by destination, its maximum physical-sector amplitude is 2.97077e-8 and squared action norm is 3.31087e-14. Its square column lies within 7.04e-15 in numerical L2 distance of the multiplier span. It is not exactly zero as a rational polynomial on the sector; these are finite-precision redundancy diagnostics, not a proof that the ray can be discarded without accounting for its effect. The diagnostic retains individual-monomial counts under explicitly unaggregated names.
