# V6: the proposed method failed its headroom gate

The privileged controlled-history representation gave **no predictive advantage over strong stateful ARX**. The two selected the same coefficients and produced forecasts equal to floating-point tolerance across 32 simulated coupled thermal systems and an independently measured TCLab record. The requested two-generation method-learning pass has **not succeeded**. Autonomous acquisition, m1→m2 transfer, and net method-learning benefit were not run because the prerequisite headroom gate failed.

V1–V5 remain regression fixtures. V5's complementary acquired-information theorem is preserved; it is not relabelled as learning a research method.

## The decisive comparison

The [protocol](../../research/v6/PROTOCOL.md) was hash-frozen before fitting. The proposed privileged procedure constructs a state from input/output history, fits an update and compiles it into a companion-state model. It receives no true state, hidden order, coefficients or poles. The conventional ARX arm can select the full declared order range using exactly the same observations and fitting/selection rules.

[The algebraic argument](../../research/v6/history_closure_theorem.md) explains the result. Cayley–Hamilton supplies finite input/output closure for noiseless LTI dynamics. Separately, for any fitted ARX coefficients, its companion realization and scalar recurrence generate the same forecasts by induction. The numerical implementations are separate; their maximum normalized forecast difference was **1.56×10⁻¹⁴**. This removes an expressivity advantage for this candidate, not for every possible research method.

The trial used independently generated networks of 3–6 hidden thermal states, two heater inputs and one measured output. White and correlated observation-noise regimes were specified in advance. Each system supplied distinct fit, selection, calibration and test trajectories. Tests used new initial conditions and pulse/sinusoidal inputs. Statistical sampling units are the 32 systems, not their thousands of within-system observations.

| Method | Mean normalized forecast RMSE | Mean block coverage | Systems meeting all quality requirements |
|---|---:|---:|---:|
| Privileged history-state compilation | 0.030958 | 84.23% | 11/32 |
| Strong stateful ARX | 0.030958 | 84.23% | 11/32 |
| Full-history retrieval | 0.490253 | 86.36% | 0/32 |
| Sparse linear-history discovery | 0.033024 | 84.09% | 12/32 |
| Current-state diagnostic | 0.204526 | 85.09% | 5/32 |

The fixed quality target was normalized RMSE≤0.15, empirical block coverage≥90%, and no validation refusal. The useful improvement over current-state prediction does not establish headroom against the strong history baseline. No positive statistical bound on an improvement is claimed: the privileged/ARX prediction equality is algebraic. The sparse model remains a competitive cheaper baseline, despite a slightly larger average RMSE.

Both privileged and ARX selected among 768 fitted candidates over the simulated systems. Their normal-equation work proxies are identical: 66,444,288. Dense companion forecasting uses a proxy of 12,345,344 arithmetic operations versus 546,304 for the scalar recurrence. These are transparent implementation work counts, not machine-instruction measurements or a universal lower bound. A sparse optimized companion implementation could match the recurrence; this would still be ordinary equivalent computation. Measured fitting times were approximately 2.33 seconds for each, with single-run timing uncertainty.

## Independent measured anchor: a separate validity failure

Dowling Lab describes real step and sinusoidal tests on its two-heater TCLab apparatus and publishes both input/output records. We froze the two 901-row files at a repository commit and retained their BSD license and hashes. [Primary experiment source](https://dowlinglab.github.io/pyomo-doe/notebooks/tclab-model/), [local provenance](../../research/v6/primary_sources.md).

Models fit the first part of the step experiment; later step segments were used for selection and interval calibration. The untouched sine experiment was the offline prediction target. All arms saw T1 and both heater inputs; T2 was not supplied. This is one apparatus and a forcing change, not cross-apparatus validation or new active experimentation.

| Method | Sine-test normalized RMSE | Block coverage |
|---|---:|---:|
| Privileged history state | 0.157650 | 3.64% |
| Stateful ARX | 0.157650 | 3.64% |
| History retrieval | 0.374406 | 18.18% |
| Sparse history | 0.146145 | 0% |
| Current state | 0.146315 | 0% |

Every method failed the joint quality requirement. The empirical intervals do not carry an unconditional coverage guarantee under serial dependence and changed forcing.

The [post-run validity diagnosis](../../research/v6/validity_diagnosis.md) identifies a concrete information gap: both training heater inputs are constant, so their centered contrast rank is zero. A change in input gain can be offset by a change in intercept without changing any training prediction, yet produce a different response to the sine input. All 55 test blocks contain an unsupported input contrast. That diagnostic uses inputs only, but was added after the frozen run and is **not** counted as successful prospective refusal. Sampling jitter and unmodelled physical dynamics are additional possible mismatches.

## Cost ledger and audit

Resources remain separate:

- **Simulation:** 320 generated trajectories, 61,440 transitions, and 61,760 scalar output records, reused identically by every comparison arm. Inputs and fitting records are also retained. These are not physical measurements saved.
- **Physical measurements newly performed:** zero. Two existing measured files were downloaded for offline prediction; no active experiment-selection advantage is asserted.
- **Repeatable computation:** every fitted candidate, sparse refit, retrieval selection, forecast and validation pass is counted in [the results](headroom_results.json). Fit/evaluation wall times and work proxies are separate fields. Failed validation remains visible. One full replay is a separate verification expense, not secretly part of every deployment's fitting cost.
- **Method acquisition/training:** not attempted. Zero new method artifacts does not mean a successful zero-cost acquisition; this gate was not reached.
- **One-off development:** source research, proof development, implementation, debugging and review are distinct from algorithm execution. The goal's recorded agent/time usage is reported separately in the conversation; it is not charged as every system's training cost.

The V5 import audit confirms **no dependency on V4's 576,000-labelled-episode calibration runs**. Those are separate historical experiments. V5's recorded 23,296 source executions comprise its transfer and retention studies; they should neither absorb the unrelated V4 episodes nor be replaced by only the final 21-read deployment figure.

**113 tests pass**, including all prior regression tests, exact forced-recurrence arithmetic, independent companion/ARX forecasts, future-output isolation, the constant-input obstruction and the V4/V5 dependency check. The [verification receipt](verification_receipt.json) reproduces the full frozen run and verifies preserved V5 source/result hashes. The full two-generation pass remains unproved.

## Gates and next decision

| Gate | Outcome |
|---|---|
| Privileged method provides headroom over strong stateful baseline | Failed for this method/family |
| Reliable physical transfer under the declared interval target | Failed |
| Autonomous acquisition of m1 | Not run |
| m1 causally reduces acquisition of m2, including disable/restore/irrelevant controls | Not run |
| Two frozen downstream method improvements | Not run |
| Positive acquisition-inclusive net method benefit | Not tested |
| Cross-apparatus or active physical experiment savings | Not tested |

The bottleneck here is **candidate-method headroom**, with a separate **validity/excitation failure**. There is no evidence yet about whether an autonomous learner could acquire or select a helpful method, because a helpful candidate has not passed this gate. A new proposal should target work the strongest existing solver actually performs inefficiently and must first pass a fresh positive-control comparison. Changing this test after seeing its outcomes would not supply that evidence.

## Reproduce

```bash
python3 -m unittest discover -s tests -v > results/v6/test_log.txt 2>&1
python3 -m experiments.v6_headroom
python3 -m experiments.v6_audit
```

Python 3.12+ and the existing NumPy dependency suffice. No new framework, model training service or library-learning package was installed. N4SID, DreamCoder, Stitch and AI Feynman were researched as relevant comparisons; none is falsely reported as an executed arm. The sparse-history implementation is a limited supplied linear dictionary, not a full SINDy implementation.
