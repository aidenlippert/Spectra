# Fixed-factor ideal repair

`ideal_repair.py` tests whether the loss from a truncated SOS factorization is
partly a representation-gauge problem. It holds every integer factor entry
fixed, constructs the exact CAR coefficient map, and solves an LP over the
number-conserving multiplier `X` and scalar `b`:

\[
  \min_{X,b,r}\;\|r\|_1-b,\qquad
  r = H-b-\sum_i L_i^\dagger L_i-(\hat N-N)X.
\]

The LP is only a proposal step. Coefficients are rationally rounded and the
result is accepted only after `experiments.marginal_symbolic.verify` rebuilds
the CAR polynomial exactly. The verifier reports the sound lower bound
`b-eta`; no sector matrix or full-sector enumeration is used.

Example bounded smoke run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  research/certificate_scaling/ideal_repair.py \
  --input results/certificate_scaling/low_rank_compression/campaign/rank_full.json \
  --output results/certificate_scaling/ideal_repair/rank_full_repaired.json \
  --receipt results/certificate_scaling/ideal_repair/rank_full_repaired_receipt.json \
  --rounding 1000000
```

The corrected rank-cap inputs under
`results/certificate_scaling/low_rank_compression/campaign_corrected/` have
verifier-compatible shapes and were replayed at ranks 8, 32, 64, and full. The
repair reduces residual l1 at every tested rank, but the verified lower bound
is not guaranteed to improve: it improves from -10.1515 to -8.7366 (rank 8),
-3.6177 to -3.4947 (rank 32), and -3.3561 to -3.3490 (rank 64), while the full
case changes from -3.33039 to -3.33048 after rational rounding. The repaired
JSON is larger because it stores the dense multiplier explicitly. Exact
comparisons are in `results/certificate_scaling/ideal_repair/comparison.json`.

An earlier uncorrected campaign artifact was malformed: for example,
`campaign/rank_2.json` contained a block with 2 factor rows and 4 words (and
other blocks with 2 rows and 24 words). It failed the verifier's required
`len(row) == len(words)` check. The current script verifies the source first
and rejects such shapes; it performs no transpose or guessed repair.
