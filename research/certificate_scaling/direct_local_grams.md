# Direct locality-only Gram proposals

`direct_local_grams.py` constructs Gram blocks from mode hyperedges appearing
in the Hamiltonian polynomial. It does not read source certificate factors or
their support. Each block contains same-charge annihilation or creation words
from subsets of one Hamiltonian edge, capped by `--max-words`; the existing
particle-number ideal and exact symbolic verifier are reused.

Bounded runs completed under `results/certificate_scaling/direct_local_grams/`:

- Square H4, `--max-words 16`: 380 local blocks, maximum block dimension 6,
  4,048 Gram scalar variables, build 0.41 s, solve 2.91 s, exact lower
  `-4.4584757993`, residual l1 `2.37e-6`, residual maximum degree 6.
- Rectangle H4, `--max-words 16`: 212 local blocks, maximum block dimension 6,
  2,200 Gram scalar variables, build 0.36 s, solve 0.73 s, exact lower
  `-4.9678181531`, residual l1 and higher-body terms are recorded in its
  receipt. The source CAR upper candidate is `-4.4758966584`, so this local
  proposal is not yet a useful tight bracket.

The square `--max-words 4` attempt was infeasible and is deliberately not
reported as an accepted certificate. This is useful evidence that locality
alone needs a larger per-edge dictionary or a better support-selection rule.
All accepted outputs were passed through exact CAR verification; no full-sector
enumeration was used.
