# Proof-carrying adaptive experiment policies

This proposal gives a finite, exact test of whether two retained scientific
mechanisms interact constructively when identifying a fresh hidden physical
parameter. It is deliberately narrower than universal law discovery: the
model class, prior, experiments, deterministic outcomes, and resource budget
are all finite and explicit.

## Finite deterministic problem

Let `M` be a finite set of hidden models, with rational prior `p(m)`. An
experiment `e` has a finite outcome map `y = f_e(m)` and rational cost `c_e`.
After history `h`, the posterior support is

`S(h) = {m : f_{e_i}(m)=y_i for every observation in h}`.

The goal is to identify one model exactly. A state is therefore the remaining
support `S`, together with the set `E` of experiments already consumed when
experiments are one-use resources. If experiments are reusable, omit `E`.

The exact optimal expected-cost recurrence is

`V(S,E) = 0` if `|S| <= 1`,

`V(S,E) = min_{e notin E} [ c_e + sum_y P(y|S) V(S_{e,y}, E union {e}) ]`

otherwise, where

`P(y|S) = sum_{m in S: f_e(m)=y} p(m) / sum_{m in S} p(m)` and
`S_{e,y} = {m in S : f_e(m)=y}`.

The recurrence is finite dynamic programming because each useful experiment
strictly refines support. A policy certificate contains, for every reachable
state, the selected experiment, every successor support, each rational branch
probability, and the claimed rational value. A verifier checks partitioning,
probability normalization, recurrence equality, and terminal singleton states.
Backward induction then proves both policy soundness and optimality; an
exhaustive lower bound is the same Bellman value, not a benchmark against a
heuristic.

For worst-case cost replace the expectation by

`W(S,E) = min_e [c_e + max_{y:P(y|S)>0} W(S_{e,y}, E union {e})]`.

This produces a separate certificate. Expected and worst-case claims must not
be conflated.

## Retained-knowledge interaction test

Define four conditions using the identical `M`, prior, experiment menu, and
costs:

* `none`: the policy receives the base model family only;
* `A`: it receives mechanism A as a fixed, reusable restriction or derived
  outcome map;
* `B`: likewise for mechanism B;
* `AB`: it receives both mechanisms and their claimed composition.

Each condition gets its own exact Bellman certificate. Report

`Delta_A = V_none - V_A`, `Delta_B = V_none - V_B`,

and the interaction residual

Define complementary savings as `J = V_A + V_B - V_none - V_AB`.

Positive interaction means `J > 0` for savings in expected cost (and it must be
tested independently for worst-case cost). The sign is not assumed. A negative
or zero result is a valid outcome. To avoid changing the task while changing
knowledge, all four runs must use exactly the same hidden models, prior,
available experiments, outcome semantics, and resource accounting; only the
retained mechanism information may differ.

The comparator must be the exact optimal posterior policy under each condition,
not a frozen baseline or a policy that is denied the learned logic by
implementation choice. If A or B changes the admissible experiment set, that
change is a separate intervention and must be reported rather than folded into
knowledge reuse.

## Rational verifier schema

A compact certificate can use records of the form:

```json
{
  "condition": "AB",
  "models": ["m0", "m1"],
  "prior": {"m0": "1/2", "m1": "1/2"},
  "experiments": {"e0": {"cost": "1", "outcomes": {"m0": "0", "m1": "1"}}},
  "states": [{
    "support": ["m0", "m1"], "used": [], "value": "1",
    "action": "e0",
    "branches": [{"outcome": "0", "probability": "1/2", "support": ["m0"]},
                 {"outcome": "1", "probability": "1/2", "support": ["m1"]}]
  }]
}
```

The verifier rejects duplicate or missing states, successor supports that do
not equal the declared outcome partition, non-rational values, probabilities
that do not sum to one, actions unavailable under `used`, and values that do
not satisfy the Bellman equality. It also checks every terminal support has
value zero. The model enumeration and experiment table are part of the
certificate hash so a policy cannot silently change its task after synthesis.

## Noisy finite-shot extension

For finite outcomes with known rational channel `Q_e(y|m)` and `n_e` shots,
replace a deterministic outcome with a count vector `z`. The exact transition
probability is the multinomial rational

`P(z|m,e) = n! / product_i z_i! * product_i Q_e(i|m)^{z_i}`.

The posterior update is rational Bayes, and the same Bellman recurrence sums
over count vectors. State support alone is no longer sufficient: retain the
posterior weight vector (or an exact canonical numerator tuple) and remaining
shot/experiment budget. With bounded shot counts this is still a finite DAG if
each action consumes at least one shot. A verifier checks the channel rows,
multinomial normalization, Bayes update, and Bellman equality exactly.

If outcomes or likelihoods are irrational, interval certificates must replace
exact rational equality and report an uncertainty gap. They cannot be silently
rounded into an exact proof.

## What this establishes and what it does not

The result proves optimality only for the declared finite model class and
experiment resources. It does test compounding reuse: AB must solve the same
fresh identification task after retaining both mechanisms, and the measured
interaction residual must be positive under the shared accounting. It does not
prove that the mechanisms are physically correct, discover mechanisms outside
the declared class, or establish a universal scientific-intelligence theorem.

The most useful failure cases are deliberate: an AB representation that gives
no savings, a mechanism that improves one posterior branch while worsening
another, and a misspecified model class where every policy is confidently
wrong. These should remain visible in the certificate and benchmark report.
