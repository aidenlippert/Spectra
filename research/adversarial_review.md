# Final adversarial review receipt

Audited `experiments/run.py` and `experiments/discovery.py` read-only. Current code hashes:

```text
experiments/run.py          a0d1dd65cf7a12c8f3f19d8035136c1b8ef941bf3a26ebf30be0373985248c17
experiments/discovery.py    3b9b7aaf938f0e90916b923f9ea683f145dcf8e8214d183cfee5b4800fec467c
experiments/pauli.py        a89313603d3a6d4589f6188df454b679204e2d80584403837dc3bb49732055ba
experiments/certificates.py f83cba5368140dd8d3992ffa354680e9c0e196da8152eb0ac76e9f01dd8b993a
```

Test receipt: `python -m unittest discover -s tests -v` ran 29 tests; all passed. The tests cover exact commutators, actual-oracle sensor agreement, three-stage support transfer, uncertainty-aware dynamics, rational certificates, independent dense algebra checks, withheld-coupling omission, and the corrected scaling cap.

## Verified current interpretation

The integrated pipeline now checks the learned model against actual dense Schrödinger trajectories from an oracle Hamiltonian, propagates the coefficient uncertainty into a conditional observable-vector bound, and requires the result to remain within that bound. Certificate construction and rational checking are separate from dense numerical validation. The CLI pipeline gate is explicit and JSON output is serializable.

The reported 18 contexts are evaluated contexts, not 18 statistically independent or unique Hamiltonians. They comprise repeated base public Hamiltonians, staged support changes, and fixed coefficient rescalings. This is same-family support/coefficient transfer within a public finite Pauli vocabulary. It is not open-ended mechanism discovery, new-family transfer, or evidence of compounding scientific intelligence; the code and report correctly keep that verdict false.

## Residual limitations (no blocking bug found)

- The derivative oracle still abstracts away physical preparation, finite-shot measurement, calibration, and apparatus cost.
- Fixed public probes mean the system does not learn an intervention or experiment-selection policy.
- The energy interval is conditional on the declared dictionary and sensor-error bound; it is not a universal physical certificate.
- The reported sequential/scratch arithmetic comparison is a partial greedy correlation cost, not end-to-end runtime or total experimental cost.
- The oracle trajectory check validates the finite declared simulator and uncertainty model; it does not establish causal closure under arbitrary new couplings.

No blocking correctness defect was found in the current bounded simulator or its 29-test receipt. The remaining limitations are scope statements, and the current report does not overclaim them.
