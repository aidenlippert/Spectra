# Adversarial requirements for a finite-domain compounding claim

This document defines the gate before results are inspected. “Compounding” means that two retained discoveries jointly reduce the cost of a withheld task more than either discovery alone, under matched information and resource budgets. It does not mean universal scientific intelligence or superiority over every possible algorithm.

## Scope and theorem boundary

The claim must name a finite task family, hidden variables, admissible observations/interventions, noise model, error target, and resource units. A constructive theorem may prove compounding for that family. The result must separately label which quantities are theorem guarantees, simulator assumptions, and empirical estimates. A finite-domain theorem is not evidence that the effect extends to unrestricted physical law discovery.

## Predeclared positive gate

Before inspecting outcomes, specify thresholds and seeds. A positive result requires:

1. Two discoveries are acquired from observations without hidden labels or answer injection.
2. Each discovery transfers to independently generated withheld instances.
3. The joint agent meets the target error/resource threshold.
4. The mixed complementarity effect is positive under a declared cost convention. For target cost (K), report (K(\varnothing),K(A),K(B),K(A,B)) and the joint gain relative to the isolated gains. Do not call a lower joint cost “compounding” if one discovery alone already gives the same benefit.
5. The result replicates over predeclared seeds or has an analytic proof covering the randomization.
6. Any claimed confidence interval or error guarantee uses only training data for adaptation and independent data for evaluation.

A null result is valid when the joint gain is within the predeclared equivalence margin, transfer fails, or acquisition cost eliminates the apparent gain.

## Four-way ablation and fair baselines

Run the same algorithmic family with mechanism state reset or retained:

| Condition | Discovery A | Discovery B | Retained state |
|---|---:|---:|---|
| neither | no | no | empty |
| A only | yes | no | A |
| B only | no | yes | B |
| both | yes | yes | A+B |

The “B only” condition must receive the same observation budget and source access as the joint condition except for A’s retained artifact. The “from scratch” baseline reruns discovery with the same code and stopping rule; the retained agent may reuse only artifacts actually acquired. Include frozen, retrieval-only, conventional fixed-feature, and an informed optimal-design baseline when the finite theorem grants it the same observations. Do not require beating an optimal baseline; a theorem of complementary reuse can hold when the optimal baseline ties the protocol.

Report acquisition, training, discovery, verification, reuse, experiment, and failure costs separately, then report cumulative cost over (N) withheld tasks. Cache lookup and artifact loading are charged. State the amortization break-even point. Equal sensor records alone do not imply equal total cost.

## Leakage and simulator boundary

Hidden mechanisms must be absent from filenames, ordering, metadata, candidate restrictions, and stopping signals. Randomize candidate order and source order. The simulator may expose an ideal oracle only if the oracle is named as an assumption; it cannot be cited as physical evidence. If finite-shot trajectories are used, sample them from the declared state preparation and measurement model, include calibration/readout noise, and preserve a separate exact-algebra checker. Do not use the hidden Hamiltonian, exact ground-state energy, or withheld outcomes in model selection, certificate construction, or threshold tuning.

## Statistical and mathematical validation

Use independent held-out instances, not merely new rows from the same generated instance. Report per-seed outcomes, confidence intervals or concentration bounds, and the number of adaptive trials. If the guarantee is analytic, show that its assumptions cover the generated test distribution. Do not turn a maximum over finite holdouts into a supremum claim.

For certificates, distinguish an exact symbolic/rational checker from dense numerical validation. A primal feasible approximation or a fitted coefficient rounded to a rational is not a certificate unless the checker verifies the claimed inequality. For dynamics, state whether the bound covers nominal model error, omitted operators, controls, and integration error.

## What counts as reuse

The retained object must be a mechanism, representation, experiment rule, or certificate procedure whose preconditions and transfer map are recorded. Reusing a support list or cached answer is retrieval; it becomes a compounding scientific result only when the retained object applies to a withheld context with lower measured acquisition or computation cost and without expanding the feature budget for free. A new coefficient vector on the same public Hamiltonian family is valid finite transfer, but must not be described as new-family mechanism discovery.

## Required report language

The final report must state: task-family scope; four ablation results; baseline access; acquisition and amortized costs; leakage controls; statistical or theorem basis; exact-check status; and the strongest unsupported extrapolation. It must explicitly distinguish a mathematically constructed complementary example from evidence about general physical science. Claims of “compounding intelligence” should be reserved for a positive mixed effect that survives these controls; otherwise report the measured finite-domain transfer or a null result.

## Review gate for the revised GF(2) mask model

The proposed model is mathematically valid as a finite-domain complementarity theorem if its assumptions are explicit: the state family is \(\rho=(I+sP_a\otimes P_b)/4\) with declared axes, local actions are the stated four operations, signs are sampled uniformly, and the response laws are linear masks \(a(x)=m_A\cdot x\), \(b(z)=m_B\cdot z\) over unseen binary control contexts. A calibration qubit observed through a binary symmetric channel may identify the masks, after which compilation selects the correct measurement. The claimed risks must be derived from the specified channel and Bayes prior, with the exact convention for \(R_0,R_A,R_B,R_{AB}\) stated before simulation. Under the supplied formulas, \(J=R_A+R_B-R_0-R_{AB}=\nu/8\) is a theorem only for that menu, prior, noise model, and finite horizon.

This positive result establishes learned-rule complementarity, not superiority over a conventional Bayesian agent with the same calibration observations. Such a baseline can infer the same posterior and tie the compiled policy. The report must separate the theorem, the algorithm comparison, and the scientific interpretation. The theorem concerns the finite response-law family; it does not imply broad scientific intelligence.

To rule out trivial answer retrieval, draw masks from a hidden distribution and evaluate on control contexts whose input bits were never seen during acquisition. Exact mask/context pairs must be absent from metadata, ordering, filenames, and prompts. An exact-answer cache baseline should fail on these contexts. A structured retrieval baseline that stores or reconstructs the mask law is allowed to tie; that demonstrates reusable structure, not unique learner superiority. Report an erasure ablation that deletes the learned mask while preserving raw calibration data, and a conventional Bayesian baseline receiving the same raw data and prior.

An independent verifier should check the response-law algebra, posterior or finite-horizon dynamic program, risk expressions, and compilation policy. Use independent seeds for withheld masks and contexts. Report analytic risk and Monte Carlo confidence intervals; do not infer the theorem from simulation alone. The valid conclusion may be: finite-domain complementary learned rules verified; no algorithmic superiority shown.
