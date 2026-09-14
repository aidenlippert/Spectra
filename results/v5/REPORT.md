# A bounded compounding result

Two pieces of information acquired from noisy measurements have **strictly complementary value** for a fresh coupled-system task. Retaining them also supports two successive measurement savings after their acquisition costs are charged. This is proved for a declared four-dimensional Gaussian instrument family and checked in executable simulations. The model grammar and experimental update rule are supplied; the result is not unrestricted scientific intelligence or a new physical law.

## What is proved

Write R∅, R_A, R_B and R_AB for optimal one-scalar prediction errors on the same fresh target, with neither, either, or both acquired sign bits retained. The proof establishes

\[
J=R_A+R_B-R_{AB}-R_\varnothing
\ge\frac{963696138679}{12369000000000}
\approx0.0779122>0.
\]

Equivalently, the two discoveries together provide more downstream error reduction than the sum of their separate benefits. The bound includes actuator uncertainty, bounded readout shift, finite-preparation tolerance and errors in the two acquisitions. It is a lower bound, not an exact value of optimal risk or a finite-sample confidence interval.

Separately, an adaptive investigator lacking the predecessor information needs at least **1.34637 expected scalar reads** to reach error 0.05 in the declared class. Acquiring the predecessor once and then solving four new targets costs **five reads total**, or 1.25 per target. This acquisition-inclusive gain occurs at two successive stages. The comparisons permit adaptive stopping and arbitrary scalar postprocessing, and grant more informative analog outputs to the lower-bound comparator.

[The full proof](../../research/v5/PROOF.md) derives the geometry, Gaussian separation, expected-read obstruction, robustness, finite relaxation and common-task complementarity. [Exact arithmetic certificates](complementarity_results.json) verify all three potentially troublesome mismatch cases. These are ordinary mathematical arguments and rational checks, not a machine-checked proof in a proof assistant. No claim of novelty relative to the entire literature is made.

## Measured transfer

The learner receives binary source readings and public contracts. Unknown modes, covariance matrices, target signs and evaluator seeds are not passed to it. A learned mode supplies a direction for nulling nuisance fluctuations in the next measurement. Every new target has a fresh component identity and coupling scale; the family grammar and shared predecessor structure remain fixed.

Across 64 independent simulated root worlds, each method faces four second-stage and sixteen third-stage targets:

| Method | Reads per world, including its acquisition | Stage 2 errors / 256 | Stage 3 errors / 1024 |
|---|---:|---:|---:|
| Cumulative acquired modes | 21 | 0 | 9 |
| Frozen after first mode | 37 | 0 | 9 |
| From scratch | 56 | 0 | 9 |
| Exact-key retrieval with scratch fallback | 56 | 0 | 9 |
| Structured replay of earlier evidence | 21 | 0 | 9 |
| Stateful Bayesian inference | 21 | 0 | 9 |
| Fixed uninformed one-read policy | 20 | 56 | 359 |

The stateful and structured-retrieval ties are substantive. The gain comes from acquiring and composing physical information, not from a demonstrated advantage over conventional learning. The two/three-read reset policies are implemented upper bounds; their optimality is not claimed. The universal cold lower bound is the weaker 1.34637 figure.

The fixed uninformed policy uses fewer reads but fails the target accuracy. Frozen and scratch are explicit retention ablations; a baseline permitted to retain the same acquired information can match the cumulative method. Exact task-key hits are zero. Full saved records and counts are in [transfer_results.json](transfer_results.json).

## Common-target retention experiment

A separate 128-world study acquires two bits before drawing fresh target signs, then gives each policy only its permitted subset of those bits. In particular, B-only does not receive a full mode vector or historical probe that reveals A. Across 2,048 fresh targets:

| Retained information | Observed policy error |
|---|---:|
| Neither | 0.375000 |
| A only | 0.232910 |
| B only | 0.250977 |
| Both | 0.002930 |

The observed policy contrast is 217/2048≈0.105957. It is consistent with the theorem but is not an estimate of optimal risks. Shared prefixes and paired randomness create dependence; no iid confidence interval is asserted. The study charges 256 calibration reads and 8,192 counterfactual target reads. See [the complete retention study](complementarity_results.json).

## Resources and validation

Every source read includes a reset, four coordinate preparations, finite relaxation, a bounded unit-vector command and a binary record. The model uses covariance scales up to about 2·10¹⁴ and relaxation **6,000,000,000,000,030 time units per read**, with actuator error at most 10⁻⁹. These are enormous physical requirements. Reset/source availability and their finite common costs are declared primitives; an actual apparatus, numerical hardware energy cost and finite-range bath implementation are not derived.

The exact OU law has a proved finite-relaxation TV bound, conditional on its coefficients and initialization. NumPy simulation does not itself certify a continuous Gaussian law. Gaussian amplitudes are unbounded; the proof gives finite expected energy and an endpoint tail bound, not a deterministic apparatus-range or pathwise guarantee.

The [computational profile](computational_accounting.json) includes fixture generation, all seven methods, both studies and serialization: 73,337,133 profiled calls, 23,296 scalar-source executions, about 61.4 CPU seconds with profiling and 80.9 MB peak traced Python allocations. Native memory, startup, earlier research and verification reruns are explicitly excluded. Actual retained rational mode coordinates use at most eight bits per numerator/denominator. The fixed-grammar learner also has a conservative polynomial bound in the finite contract description length; no general physical-solver complexity result is claimed.

**106 tests pass.** Tests include exact mode geometry and residual Gram identities, source and error arithmetic, refusal before invalid experiments, evidence corruption, restricted B-only access, and an intentionally undeclared control mismatch that destroys the promised accuracy. A stored contract does not automatically detect arbitrary real-world misspecification.

The [verification receipt](verification_receipt.json) records replay of 14,848 saved transfer source reads, 8,960 target decisions, 6,400 reconstructed learned artifacts and the full retention experiment. Replay establishes consistency and reproducibility; the independent mathematical argument establishes the probability bounds.

## What this changes

V4's trained-policy experiment was a null result. It exposed a missing joint-information requirement. V5 supplies a concrete example in which acquired physical information changes which experiment is useful, transfers twice with net cost savings, and has strictly complementary downstream value. This meets a bounded operational definition of compounding. It does not show that the system learned the representation grammar, discovered a new interaction vocabulary, transfers to general materials, or solves the ten keystone endpoints.

[The requirement audit](../../research/v5/goal_audit.md) separates the satisfied bounded milestone from those remaining ambitions. The immediate next scientific escalation would remove the supplied recursive mode grammar while preserving the same certificates and fair comparisons.

## Reproduce

From the workspace with Python 3.12+ and NumPy:

```bash
python3 -m unittest discover -s tests -v > results/v5/test_log.txt 2>&1
python3 -m experiments.v5_run
python3 -m experiments.v5_complementarity
python3 -m experiments.v5_account
python3 -m experiments.v5_verify
```

No physical experiment, remote publication, commit, paid API call or external message was performed as part of this experiment.
