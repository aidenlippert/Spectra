# Boundary for a cumulative scientific-capability claim

A growing library of examples is not evidence that a research system has
acquired a reusable scientific operation. The finite claim that can be tested
is conditional value of information on a sequence of fresh tasks.

Let `D_t` be the retained discoveries after stage `t`, and let `T_{t+1}` be a
fresh hidden-model task drawn from a declared distribution. Fix the same
resource budget, model class, prior, likelihood table, and admissible actions
for all systems. Let `R(A,T)` be exact Bayes risk plus a declared resource
penalty for algorithm `A` on task `T`. Define the conditional value

`V(D_t) = E[R(A_base,T_{t+1}) - R(A_{D_t},T_{t+1}) | D_t]`.

Discovery `d_{t+1}` compounds capability only when the held-out increment is
positive after charging acquisition cost:

`Delta_t = V(D_t union {d_{t+1}}) - V(D_t) - cost(d_{t+1}) > 0`.

The test must use a fresh task family and a predeclared split. Replaying the
examples used to create `d` measures memorization or amortized lookup.

## Conventional tie is compatible with a positive cumulative claim

Suppose a conventional algorithm `C` has the same risk as the learned system
on every task in the current finite benchmark. That establishes only

`R(A_{D_t},T)=R(C,T)` on that benchmark.

It does not imply `Delta_t <= 0`. A reusable representation can tie an
emulator's outputs while reducing the number of interventions, proof checks,
or experiments required on a new family. To make the stronger claim honestly,
report both predictive regret against `C` and conditional value under a second
resource axis (experiment count, certified computation, or transfer cost).

If the metric is only terminal prediction risk and `C` is itself Bayes-optimal
for the same task distribution and budget, then no positive risk improvement is
possible: Bayes optimality gives `R(A,T) >= R(C,T)` in expectation. A claimed
gain must therefore concern a declared resource or a shifted fresh task
distribution, and that shift must be fixed before evaluation.

## Sequential reuse test

Use stages `1,...,q`. At each stage, train or discover on family `F_t`, freeze
the retained artifact, then evaluate on unseen family `F_{t+1}`. Compare:

* a stateless baseline retrained from scratch;
* a stateful conventional baseline with the same storage and compute budget;
* the cumulative system reusing `D_t`.

Charge discovery, labeling, teacher-planning, storage, and verification costs
separately from reuse cost. Require `Delta_t>0` for at least two successive
stages with confidence intervals or exact finite enumeration. A single positive
transfer can be a lucky task family; repeated positive increments are evidence
of compounding reuse within the declared domain.

## What cannot be proved from finite experiments

No finite benchmark proves an unbounded theorem that the system will keep
improving. A mechanism may stop transferring, and an adversarial fresh family
can make every retained representation unhelpful. The strongest honest result
is a finite theorem conditional on the enumerated model family: exact Bayes
risk, exact resource accounting, and a positive held-out conditional value.

Thus “the system gets better at science” should be stated as a measured,
resource-accounted sequence of positive conditional values, with explicit
failure and abstention cases—not as superiority over a conventional emulator
that already ties it on the current prediction benchmark.
