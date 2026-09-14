# Matched CAR-SOS comparison

Run with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.certificate_scaling.operator_pricing_comparison
```

The receipt is [receipt.json](/Users/aidenlippert/Documents/Spectra/results/certificate_scaling/operator_pricing_comparison/receipt.json).
This is a 4-mode, 2-particle Hubbard Hamiltonian with seven mixed CAR Gram
blocks. It compares the full dictionary with fixed, seeded-random, and greedy
adaptive subsets at matched budgets 1, 2, and 3. Every greedy candidate score
is recorded, including infeasible statuses and elapsed pricing rounds.

The symmetric variational physical upper estimate is approximately
`-0.1403124221`. All feasible strategies in the corrected run exported and
passed exact replay. Certified gaps to that upper estimate ranged from about
`3.81e-7` (fixed, budget 3) to `2.26e-6` (adaptive, budget 3). Budget-1 and
budget-2 fixed selections were numerically infeasible; this is recorded as
solver status and is not treated as an exact cone-separation result.

An earlier receipt incorrectly reported export rejection because it attempted
to convert a rational string to `float` while serializing the comparison
report. That reporting bug is fixed; it was not a verifier failure.
