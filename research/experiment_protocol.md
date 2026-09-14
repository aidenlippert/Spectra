# Proposed stronger protocol for intervention-stable discovery

This is the initial design proposal, not a claim that the delivered baseline implements every requirement below. The actual bounded protocol, fixed thresholds and deviations are recorded in [decision_gates.md](decision_gates.md) and [the result report](../results/REPORT.md). In particular, the delivered baseline uses three qubits, an ideal derivative oracle and fixed public probes; it does not implement shot-based intervention selection or unknown-vocabulary learning.

## Scope and claim boundary

This is a bounded benchmark for support identification, intervention-stable reduced modelling, sequential discovery, and energy certification in small spin systems. The generic Pauli basis is known in advance; therefore discovery means identifying hidden supports/terms (and reusable transformation rules), not open-ended discovery of physical law. Results are conditional on the finite system size, Hamiltonian family, time horizon, noise model, and intervention distribution.

No finite holdout certifies a supremum over all states or interventions. Report empirical maxima on a separately generated test set and call them test-set estimates. A claimed uniform bound requires an analytic argument over a declared domain.

## Environment and staged task

Use four or five qubits and a public generic Pauli-string candidate set (weight at most three). Generate sparse Hamiltonians from a hidden support and rational coefficients. Expose only preparation records, control pulses, and noisy measurement outcomes. Do not expose hidden support, term names, residuals computed with hidden terms, or the random seed.

Run three sequential stages. Stage 1 contains the public baseline support. Stage 2 adds a hidden noncommuting or longer-range term and changes the intervention distribution. Stage 3 adds a second hidden term in a new held-out context. The learner may retain operators, certificates, assumptions, and counterexamples from earlier stages. Evaluate each update on new initial states, pulses, and Hamiltonians generated from independent seeds.

Interventions must include local rotations, noncommuting quenches, pulses on previously untouched qubits, and couplings between subsystems. Include a deliberately invalid compact representation; the correct output is a failure report with a witness intervention, not forced model expansion.

## Representation task

For retained observables (O_i), fit

\[
 \dot z=A(u)z+r(\rho,u),\qquad z_i=\operatorname{Tr}(\rho O_i).
\]

Select support and memory variables using only training observations. Penalize observable count, intervention count, and estimation cost. A representation is reusable only if its operator/support rule transfers to an independently generated Hamiltonian or coupling, rather than merely retrieving the same example. The report must include training residuals, test-set residual distributions, worst tested residual, and all domain restrictions.

## Energy certification

Construct an upper bound from a supplied trial state, with its preparation and evaluation cost. Construct a lower bound from exact Pauli algebra: for example, use a rational anticommuting clique (\{P_j\}) with (P_jP_k=-P_kP_j), so

\[
 \left\|\sum_j a_jP_j\right\|\leq\sqrt{\sum_j a_j^2}.
\]

Combine such norm bounds with an explicitly stated decomposition to obtain a lower bound. A numerical SDP/primal feasible point is not itself a lower-bound certificate; every reported bound must be independently checked by exact rational or symbolic algebra, or by a separately verified dual certificate. Exact diagonalization is validation only and must not construct the bound.

Report (E_-\leq E_0\leq E_+), interval width, certificate size, arithmetic-check status, and the cost of discovering and verifying the certificate. Never silently substitute energy density for absolute energy.

## Matched baselines and accounting

Use the same data, intervention budget, numerical precision, and stopping rules for:

1. **Frozen:** stage-1 representation and certificate applied unchanged.
2. **Retrieval-only:** may retrieve prior entries but cannot invent supports, interventions, or certificates.
3. **From-scratch:** reruns discovery independently at every stage.
4. **Conventional:** fixed Pauli features, least-squares dynamics, and fixed certificate search without adaptive discovery.

Report separate (C_{train},C_{discover},C_{certify},C_{reuse},C_{experiment}), plus cumulative and marginal costs. Failed experiments and rejected hypotheses count. Cached artifacts are not free: show both amortized and full-cost curves.

## Pre-registered success and null criteria

Before inspecting results, specify thresholds for support recovery, held-out predictive error, certificate interval width, transfer cost reduction, and replication across seeds. A positive sequential result requires: (a) first hidden support is inferred without label leakage; (b) second-stage discovery improves an independent held-out context; (c) the second discovery reduces cumulative or marginal cost relative to matched from-scratch and retrieval-only baselines; (d) exact certificate checks pass; and (e) the learned rule transfers beyond the training Hamiltonian.

The result is null if the apparent gain is explained by caching, larger feature sets, more measurements, favorable seeds, or a certificate verifier that uses the exact ground-state energy. A null result is scientifically valid and should identify whether the bottleneck was support identification, intervention coverage, transfer, or certification.

## Anti-leakage and falsification checks

Use independent seeds for train, discovery, and final evaluation. Keep hidden supports out of filenames, metadata, and candidate ordering. Shuffle candidate-string order. Include adversarial interventions selected after fitting. Compare against a cache-only implementation with identical serialization and lookup costs. Require transfer to new initial states and at least one new coupling. Log every rejected hypothesis and every failed certificate.

The experiment may establish finite-sample transfer in the declared benchmark. It cannot establish universal causal closure, unrestricted law discovery, or a supremum guarantee outside the tested and analytically covered domain.
