# V12: small-cover preflight closed

The one-sweep, common-dyadic overlap construction passes the unchanged V8
certificate checker at degree 11 on the existing n=3 XXZ, gamma=2, T=1/2,
epsilon=1/1000 case. The V11 Frobenius-guarded baseline requires degree 12.
This does not pass the declared complete-cost headroom gate.

The [prospective protocol](../../research/v12/PREFLIGHT_PROTOCOL.md) required at
least 5% improvement in both aggregate and median paired complete cost, with
18/21 paired wins. The complete calculation includes model/representation
construction, rejected orders and ordinary bounds, cover construction, fresh
checking and exact endpoint evaluation. One initial call per arm was recorded
separately, followed by 21 randomized three-arm blocks.

| Arm | Total seconds, 21 calls | Median paired cost / baseline | Wins / 21 |
|---|---:|---:|---:|
| V11 Frobenius guarded | 0.220736 | 1.000000 | — |
| V12 small cover | 0.207739 | 0.993156 | 12 |
| V12 disabled hook | 0.217430 | 1.005388 | 9 |

The aggregate improves about 5.9%, but median improvement is only about 0.7%
and paired wins are insufficient. This is an inconclusive small timing effect,
not a robust algorithmic advantage. No broader benchmark or learner follows.
The implementation is preserved as a conventional candidate and closed in its
current form; no further threshold tuning is scheduled.

The construction deduplicates the union of two existing ordinary partitions,
performs one untrusted floating coordinate sweep and exports exact split
coefficients using common denominator 2^20. Integer allocation reconstructs
each Pauli coefficient exactly. Fourteen groups / 60 incidences suffice. The
checker rederives the residual with the correct negative tail sign, checks
reconstruction, anticommutation, rational norm upper bounds, horizon and total
error. Floating optimization quality cannot bypass these checks.

The earlier 61-group check-only diagnostic is historical and is not a lower
cost bound for this 14-group construction. The actual small-cover check-only
diagnostic saved about 0.17 ms; complete-cost evidence above governs disposition.

[Raw paired results](preflight.json) · [exact archive](certificate_archive.json)
· [verification receipt](verification_receipt.json).

V1–V11 remain preserved. V4 calibration remains outside V5's dependency ledger.
No acquired method, held-out evaluation, physical validation, or compounding
scientific capability is claimed. The user's subsequent mission correction
makes constructive prediction, synthesis and invention the central objective;
this solver optimization is supporting historical evidence.
