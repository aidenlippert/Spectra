# Constructive reachable preparation and operation policies

This note isolates a bounded result that can connect a desired physical behavior to a
preparation-and-operation policy. It is deliberately a promise-domain result; it does
not claim a universal compiler for matter.

## Model and specification

Fix a finite horizon `T`, a compact state set `X` (or an explicitly declared finite
truncation), compact preparation/operation input sets `U_p,U_o`, and a compact model
uncertainty set `Theta`. For each `theta in Theta`, the physical state obeys

```text
x_dot = f_theta(x,u),    y = h_theta(x,u),
```

with measurable piecewise-constant controls and a reset/preparation transition system
for feedstocks, fabrication actions, and assembly. A policy is a causal map from
observations and past actions to the next admissible action. Its cost includes all
preparation and operation resources and is bounded by `B`. The requirement `phi` is a
finite-horizon trace predicate (including terminal tolerances, safety constraints, and
resource limits). Uncertainty is adversarial for robust claims; stochastic noise can be
included by replacing sets with confidence sets and requiring the stated probability.

Sufficient assumptions for a terminating constructive procedure are:

1. `X`, action alphabets, and `Theta` are compact and effectively represented by rational
   boxes/polyhedra; the reset graph is finite and has a computable cost for every edge.
2. `f_theta` and `h_theta` are computable and uniformly Lipschitz in `(x,theta)` on each
   mode, with known constants `L_f,L_h`, and controls have bounded amplitude and a
   declared switching/grid resolution. Solutions are unique on `[0,T]`.
3. A sound one-step enclosure routine returns `Post^+` and `Post^-` for every state,
   parameter, and control cell, with Hausdorff over- and under-approximation error at
   most `e_step`; the error bound is independently checkable (interval Taylor models,
   validated ODE integration, or a certified simulator).
4. The physical interface declares all variables through which modules couple (shared
   resources, timing, memory, baths, and conserved quantities). The enclosure includes
   their cross-module uncertainty.

The Lipschitz requirement is a sufficient condition, not a physical law. Discontinuous
switching, hidden memory, and model-class misspecification must be represented as extra
state/uncertainty or the result is out of scope.

## Inner and outer reachable sets

Partition the bounded policy space into finitely many cells `C` (discrete reset paths,
control words, and parameter boxes). Propagate each cell through the validated one-step
routine. Define

```text
R^-(k+1) = union_C Post^-(R^-(k), C),
R^+(k+1) = union_C Post^+(R^+(k), C).
```

The same construction applies backward from a target set, yielding inner and outer
backward reachable sets `BR^-` and `BR^+`. Soundness is inductive:

```text
R^-_k subseteq R_k(theta,pi) subseteq R^+_k
```

for every model and policy represented by the relevant cell. A cell whose inner trace
set lies in `phi` gives a constructive policy: execute its reset sequence and control
word, then use the associated observation feedback guards. A cell whose outer trace set
is disjoint from `phi` is an impossibility certificate for every policy in that cell.
The global certificate is obtained by covering the declared budget-bounded policy
space: if every cell is outer-disjoint, no admissible policy satisfies `phi`.

For a requirement expressed by a terminal loss `J <= q` and safety functions
`g_j <= 0`, the inner test can use conservative bounds
`J^+(C) <= q` and `g_j^+(C) <= 0`; the outer impossibility test uses
`J^-(C) > q` or `g_j^-(C) > 0`. Refinement splits unresolved cells and shrinks the
integration step. A valid positive result therefore consists of the policy, a trace
certificate, and a reproducible numerical error budget; a valid negative result consists
of the covering and the cellwise separating bounds.

## Robust margin, abstraction refinement, and nonattained infima

Let `Phi_gamma` denote policies whose worst-case requirement margin is at least `gamma`:
`J <= q-gamma` and `g_j <= -gamma`. If the total propagated and abstraction error is
`e < gamma`, then an inner certificate for `Phi_gamma` realizes the original `phi`.
Conversely, if every admissible policy violates `phi` by at least `gamma`, outer bounds
with error below `gamma` certify impossibility. Instances within the unresolved band
must return `unknown`; this is essential for soundness.

An abstraction-refinement loop is therefore:

```text
cover budget-bounded policy cells;
propagate inner/outer trace bounds;
accept an inner-feasible cell or discard an outer-infeasible cell;
otherwise split the largest error/cost cell and repeat.
```

This is also an approximate-simulation interface. If an abstract transition system is
related to the physical one by an alternating `epsilon`-simulation, the abstract policy
can be refined to a physical policy while increasing output/trace error by at most the
simulation budget. Pola, Girard, and Tabuada construct finite symbolic models for
incrementally globally asymptotically stable systems with selectable approximation
precision; Pola and Tabuada extend the idea to disturbances. These hypotheses are the
reason to state stability and interface assumptions rather than silently apply a grid
to arbitrary dynamics ([Pola–Girard–Tabuada 2007](https://arxiv.org/abs/0706.0246),
[Pola–Tabuada 2007](https://arxiv.org/abs/0707.4205)).

The resource optimum need not be attained. Define
`V_B = inf{ C(pi): pi satisfies phi }`. For every `eta > 0`, run the procedure on the
budget `V_B + eta` (or on a computable descending sequence of candidate budgets). If a
policy with margin `gamma` exists at some candidate budget, the procedure eventually
returns one when refinement error falls below `gamma`. If all budgets below a claimed
threshold are outer-infeasible, that threshold is a sound lower bound. At the exact
infimum, return `unknown` unless compactness plus closed dynamics/cost and a closed
specification establish attainment. Thus the theorem promises arbitrarily near-optimal
constructive policies, not an invented minimizer.

## Complexity and what is actually guaranteed

Let `N_x` be the number of state boxes, `N_theta` uncertainty boxes, `N_u` control
cells, `N_m` reset paths, and `K=T/dt` validated steps. A direct dynamic-programming
implementation costs at most
`O(K N_x N_theta N_u N_m)` enclosure calls, times the cost of one certified integrator;
memory is `O(N_x N_theta N_m)` per layer. In dimension `d`, a uniform box mesh with
resolution `rho` has `N_x = O(rho^-d)` (and analogous factors for uncertainty and
controls), so the curse of dimensionality is explicit. Branch-and-bound may do better
on benign instances but has no general polynomial guarantee. The output is complete
only on the separated promise domain (positive or negative margin at least `gamma`),
and soundness depends on the physical model class and validated enclosures.

Bounded delta-decision gives a closely related logical formulation: for computable
bounded ODE/hybrid formulas, a procedure can return the formula or a delta-strengthened
negation, with unresolved boundary behavior made explicit ([Gao, Kong, Chen, Clarke
2014](https://arxiv.org/abs/1404.7171)). It does not provide universal decidability for
unbounded, discontinuous, or misspecified physical systems.

## Failure of weakened assumptions

* Without compactness or a finite action/budget bound, a finite cover need not exist;
  an optimizer can chase a sequence of cheaper policies with no minimizer.
* Without known uniform Lipschitz/modulus bounds, a cell's corner evaluations do not
  enclose its interior. A narrow fast excursion can invalidate both feasibility and
  impossibility claims.
* Average model accuracy is insufficient: an optimizer can select the small region where
  the model's error is favorable. A certificate must cover the adaptively selected
  policy and uncertainty set.
* Independent local certificates do not compose when a shared bath, memory kernel,
  timing variable, or conserved quantity is omitted from the interface; the coupled
  trajectory can leave the claimed reachable set.
* Approximate simulation without alternating/control-interface conditions can transfer
  trajectories but not the ability to choose matching controls, so it cannot justify a
  synthesized physical policy.
* Exact boundary decisions are generally not promised. Requiring an answer when the
  optimum lies at the tolerance boundary turns a sound refinement method into an
  unsound numerical assertion.

The resulting research target is concrete: a reusable preparation-and-operation policy
with a positive certificate or a budget-relative separating witness, plus an explicit
`unknown` outcome whenever model, margin, or abstraction assumptions do not resolve the
physical question.

## Primary sources

* G. Pola, A. Girard, P. Tabuada, “Approximately bisimilar symbolic models for
  nonlinear control systems,” arXiv:0706.0246.
* G. Pola, P. Tabuada, “Symbolic Models for Nonlinear Control Systems: Alternating
  Approximate Bisimulations,” arXiv:0707.4205.
* S. Gao, S. Kong, W. Chen, E. Clarke, “Delta-Complete Analysis for Bounded
  Reachability of Hybrid Systems,” arXiv:1404.7171.
* H. Yin, A. Packard, M. Arcak, P. Seiler, “Finite Horizon Backward Reachability
  Analysis and Control Synthesis for Uncertain Nonlinear Systems,” arXiv:1810.00313.
