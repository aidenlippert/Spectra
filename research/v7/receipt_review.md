# V7 receipt review

Scope: read-only audit of `results/v7/REPORT.md` against the saved V7 JSON/log artifacts and the checker implementation. No code, report, or result artifact was changed.

## Findings

I found no unsupported claim within the requested scope.

- `prediction_requests.json` directly supports 897 training rows, order 4, `future_response_outputs_read: false`, an exact verified witness, and prediction disagreement `2/5` °C at Q1=75. Its `eta_float` is 0.6820845904858444 °C. The report correctly describes this as a model-class/residual-enclosure result and explicitly disclaims a hardware noise bound or future-process validation.
- The checker implementation recomputes residual coefficients from the submitted rational polynomial and generator, enforces requested total time, checks record coverage, pairwise anticommutation, exact rational outward square-root bounds, and claimed-bound equality. The frozen quadratic checker additionally checks square products and full term coverage; its sum form checks finite block coverage. The implementation and `residual_math.md` support the report's conditional residual-certificate wording. The report does not promote the checker to a universal theorem.
- Row counts and feasibility counts match the artifacts: 144 initial rows; BFS 27/36; adaptive Taylor 26/36; pruning data has 72 rows representing two 36-case policies, with 26 shared successes; `strong_baseline.json` has 36 rows and 26 certified. The report's 2.36 s versus 11.49 s statement is appropriately qualified as separate development runs with unequal success counts.
- The complete-cost artifacts contain 36 cases × 7 paired repeats × 2 arms. For the order-saving XXZ case, baseline order 12 versus quadratic order 11 is explicit; full square checks use 435 square-pair tests and the partial block uses 190. The report correctly describes the observed timing as a complete-cost comparison and as a failed efficiency gate, without confidence intervals, generalization, or proof-compounding claims.
- `triple_block_diagnostic.json` contains 15 repeats per arm, with 105 square-pair tests and order 11 for the adaptive block. The report correctly says the small median difference is observed in 10/15 repeats and remains unresolved rather than headroom.
- `certificate_archive.json` has 125 entries, and `verification_receipt.json` reports 125 reparsed accepted certificates, 10 replayed strong-baseline refusals, the measured witness replay, and 22 historical hash checks. It also explicitly records no acquisition, no physical validation, and `compounding_capability: not_established`.

## Boundary note

The runtime ratios are repeated measurements on the declared development cases. They should continue to be read as descriptive artifact evidence, not as a statistical confidence claim or a proof that any mechanism compounds across cases. The current report already states that limitation and keeps the goal active; no correction was needed.
