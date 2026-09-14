# Final v2 adversarial audit

Audited `experiments/v2_run.py`, `v2_laws.py`, `v2_physics.py`, `v2_policy.py`, and `tests/test_v2_integration.py`. The bounded integration test command

```text
python -m unittest tests.test_v2_integration -v
```

passes all 6 tests. The runner is materially honest about its scope: it uses hidden masks only in the `World` and evaluator, exposes calibration/reference callbacks to acquisition, verifies exact GF(2) certificates and the finite-horizon policy tree independently, and reports that conventional structured inference ties.

## Mathematical validity

The response-law theorem is valid only for the fixed visible two-axis/four-action family, visibility (1/2), BSC calibration model, declared uniform model prior, and finite horizon. `v2_policy.py` enumerates all eight ((a,b,s)) models and verifies every action branch. The reported one-shot complementarity (J=1/16) is therefore a model-conditional identity, not a general scientific result.

The B-channel acquisition is correctly conditional on a successful A mask: it uses the A-derived compiled axis and its postprocessing channel. The runner accounts for the union-bound failure probability and labels the resulting AB risk as unconditional. A failed A acquisition can invalidate B, so this conditional dependency must remain explicit in any report.

The fixed-horizon budget is also correctly labelled. The `none,h=4` policy certificate supplies a reset target-only risk reference; the `5 * trials` comparison is a committed fixed-horizon lower-budget statement, not an expected-stopping-time lower bound and not a globally learning baseline.

## Verified anti-leakage and baselines

`acquire` receives only sensor callbacks and never a `World`; hidden masks are not passed to model selection. Held-out contexts are outside the training basis. The exact lookup baseline has zero hits in the integration test. The conventional exhaustive GF(2) maximum-likelihood fit receives the same raw records and ties the learned masks, as it should. This prevents an algorithm-superiority claim while preserving the learned-rule complementarity theorem.

## Remaining limitations

1. The fixed `contexts(d)` menu is a visible finite basis. Generalization to unseen contexts follows from the assumed linear GF(2) law; it is not discovery of the law family.
2. The exact lookup and structured-retrieval ablations are useful, but structured retrieval intentionally ties the compiled learner. Any positive conclusion must say “complementary reusable rule verified,” not “learner beats conventional inference.”
3. Costs include acquisition shots, algebraic operations, verification, storage, and simulated target tasks, but do not model laboratory apparatus, physical calibration time, or learner implementation energy. The policy certificate-generation wall time is reported separately from task costs.
4. Monte Carlo error intervals are conditional on the fixed acquired masks; the exhaustive held-out-context risk is the exact finite conditional quantity. Neither establishes performance outside the declared model family.
5. `World` is a finite simulator with idealized Born/readout sampling. `physical_experiments_performed: false` is correct.

No blocking mathematical or hidden-truth access defect was found in the current runner/tests. The correct claim is a verified finite-domain complementary acquisition theorem with matched conventional tie, not broad scientific compounding or algorithmic superiority.
