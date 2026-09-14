# Operational definitions for compounding scientific capability

This note fixes what the program can claim from sequential experiments. The current PROGRAM correctly treats its first instrument as a bounded baseline: it has not shown open-ended discovery, learned intervention choice, or compounding intelligence.

## A finite-domain protocol

Let a domain be a distribution D over publicly describable problem families x, with finite resource budget B, loss l, and an admissible certificate predicate C. An instance contains the physical description and public measurement interface, but not hidden simulator parameters or evaluation labels. Split instances into acquisition A, held-out evaluation E, and (for sequential tests) a later family E2. Fix splits and stopping rules before scoring.

For learner L, report total cost K = acquisition + search + fit + verification + reuse. Report coverage-weighted loss R = E[l(L(x), y)], with confidence intervals from independent held-out instances or a declared paired design. Independent means independent evaluation families, not repeated coefficient rescalings of one simulator.

## Four distinct claims

**Cumulative performance improvement.** A sequence L0, ..., Lt improves on a declared domain if some later stage has risk at least Delta lower or cost at least an eta fraction lower than an earlier stage, with the other quantity within a preregistered tolerance, paired confidence excluding the null, and no increased unsupported-claim or certificate-failure rate. Acquisition cost is charged to the stage that incurred it. Lower search count alone is insufficient.

**Superadditive interaction of learned facts.** Let F and G be retained facts or operations. Let V(S) be a preregistered value score for fact set S, such as negative total cost at fixed coverage and error. Their interaction is I(F,G) = V({F,G}) - V({F}) - V({G}) + V(empty). A positive result requires I >= gamma on a fresh family E2, with matched representations and equal public acquisition budgets. Bootstrap or paired-randomization intervals can test I > 0. If facts overlap, use an encoding or orthogonalization fixed before evaluation; otherwise the decomposition is arbitrary.

**Amortized cost advantage.** For tasks x1 through xn, compare average cost (acquisition + sum solve costs) / n. L has an amortized advantage over baseline B only if, after a preregistered break-even n0, its average cost is lower by eta at matched coverage and error. Report acquisition, storage, verification, forgetting, and break-even costs separately.

**Unrestricted scientific intelligence.** This is not an identifiable finite-domain property. It quantifies over every physical problem, every conventional algorithm, and every retrieval system while allowing the learner to change its own representation and experiment class. No finite test establishes it. The defensible substitute is a finite theorem: after bounded acquisition, the system discovers a reusable operation q with a machine-checkable validity domain, transfers q to a fresh family, and improves against declared comparator classes.

## Why universal superiority is ill-defined

“Better than every conventional algorithm” is not a meaningful unrestricted comparison. Any learner L can be wrapped as an algorithm A_L that simulates its update rule and memory. Thus L cannot beat every algorithm unless “conventional” is restricted by explicit resources or information. The same applies to “every retrieval system”: a retrieval system can store L's policy and acquired library or simulate its computation.

Use a closed, reproducible comparator class:

1. Scratch: rediscover or refit from current public data with no retained facts.
2. Frozen: fixed representation and procedure selected before the evaluation family.
3. Retrieval/refit: access to prior artifacts but no new discovery operation; charge retrieval, adaptation, and verification.
4. Strong conventional: a specified full model, optimizer, active-design method, or symbolic solver with identical inputs, measurement budget, precision, and checker.

The learner should beat at least one strong baseline on a fresh family while preserving coverage and certificate validity. Stronger claims require beating every member of this declared class, never an undefined universal class.

## The finite transfer gate

**Theorem.** Suppose q is selected using acquisition data A only, a checker proves C(q,x) on every x in a declared promise set P, and E2 is sampled independently from D restricted to P. If the checker accepts and paired evaluation establishes Rq <= RB - Delta with confidence 1 - alpha, then the transfer claim is valid for D restricted to P, the stated resources, and the measured loss. It does not imply unrestricted intelligence.

**Proof.** Selection is independent of E2, so the confidence statement is a held-out estimate under the declared distribution. The checker establishes only C on P. Extending either result outside P, beyond the budget, or to another distribution requires another theorem. QED.

## Counterexamples and statistical failure modes

**Label leakage.** Hidden Hamiltonian coefficients, simulator seeds, or evaluation labels entering representation selection let a system encode the answer key. The held-out estimate is invalid because E2 is not independent of selection.

**Cache masquerading as discovery.** A dictionary containing every public Pauli support can reduce later greedy search. This is retrieval/refit reuse, not a discovered operation, unless the system learns a rule that constructs or selects representations outside the supplied dictionary and pays acquisition cost.

**Noncomposing facts.** Two facts can each help alone while their joint validity domains conflict, producing I(F,G) <= 0. Sequential success therefore does not imply compounding; joint evaluation is necessary.

No-free-lunch results reinforce the boundary: Wolpert and Macready show that average superiority over all finite objective functions is unavailable without assumptions on the problem distribution ([IEEE paper](https://doi.org/10.1109/4235.585893)). Every positive claim must name the structure enabling transfer and test whether the operation survives a new coupling, intervention, and comparator.
