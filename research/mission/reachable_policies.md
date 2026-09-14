# Reachable preparation-and-operation policies (corrected)

This note gives a deliberately bounded theorem for robust, finite-horizon physical
design. The earlier draft is preserved as
`research/mission_reachable_policies.invalidated.md`; its inner/outer quantifiers were
incorrect.

## Explicit bounded model

Fix a finite mode/reset graph, a finite horizon `T`, a fixed initial state `x0`, and a
finite-dimensional state `x` (or declared truncation). A *policy candidate* is an
open-loop finite word `p=(m_1,u_1,...,m_K,u_K)` with fixed switching times, where each
mode/action lies in a compact rational box and its preparation plus operation cost is
at most `B`. Let `theta` range over a compact rational-box uncertainty set `Theta`.
For every candidate and model,

```text
xdot = f_theta,m(x,u),     y = h_theta,m(x,u),
```

has a unique trajectory on `[0,T]`. Assume computable uniform moduli of continuity
(uniform Lipschitz constants suffice), bounded state/control domains, and a validated
integrator that returns an enclosure `E(p,Theta)` containing the complete output trace
for *every* `theta in Theta`, with certified Hausdorff error `e`. The specification is a
closed finite-horizon trace set `Phi`; for a margin `gamma>0`, define `Phi_gamma` by
shrinking every safety/terminal tolerance by `gamma`. Persistent uncertainty is one
single `theta` throughout a trajectory, not independently resampled at each step.

These assumptions are sufficient conditions for the theorem, not claims about arbitrary
matter. Feedback can be handled only after adding a finite observation/history state
and proving that its policy class is finite and its interface is sound; it is excluded
from the complexity claim below.

## Sound constructive and impossible outcomes

For a concrete candidate `p`, the robust feasibility condition is

```text
forall theta in Theta: trace(p,theta) in Phi.
```

Because `E(p,Theta)` is an outer enclosure, the positive certificate is the inclusion
`E(p,Theta) subseteq Phi`. This inclusion certifies one executable policy against all
models in the declared uncertainty set. A numerical center or an existential inner
set is not a positive robust certificate.

For a policy-space cell `C` (a box of action amplitudes/times and one fixed discrete
mode word), compute an outer trace enclosure `E(C,Theta)` covering
`{trace(p,theta): p in C, theta in Theta}`. If
`E(C,Theta) intersect Phi = empty`, then every policy in `C` fails for at least one
model (indeed, no trace in the enclosure satisfies the requirement), so `C` can be
discarded. If all cells covering the budget-bounded candidate space are discarded, the
result is a sound impossibility certificate for all candidates in that space.

The two tests are therefore:

```text
realization:       E({p},Theta) subseteq Phi;
cell impossibility: E(C,Theta) intersect Phi = empty for every budget cell C.
```

The second test is a universal cover argument; it must not be confused with a backward
reachable union, which is only existential over policies. A reachable union can support
existence search, but by itself cannot certify a fixed policy or robust feasibility.

## Refinement and margins

Refine an unresolved policy cell or integration interval. If the specification has
robust margin `gamma` and the total certified enclosure error is `< gamma`, then the
outer inclusion test for `Phi_gamma` proves realization in `Phi`. If every policy
violates `Phi` by margin `gamma`, the cell-disjointness test eventually proves
impossibility once its error is below `gamma`. Boundary cases return `unknown`.

Approximate simulation can replace direct enclosures only with an explicit alternating
simulation/control interface: every abstract action must be implementable by a concrete
action and the output discrepancy must be charged to the specification margin. Pola,
Girard, and Tabuada construct finite symbolic models with selectable precision for
incrementally globally asymptotically stable systems; Pola and Tabuada treat
disturbances ([primary source](https://arxiv.org/abs/0706.0246),
[alternating-simulation source](https://arxiv.org/abs/0707.4205)). Stability and the
interface hypothesis cannot be omitted for arbitrary physical dynamics.

## Nonattained resource infima

Let `V = inf{cost(p): forall theta, trace(p,theta) in Phi}` over the declared compact
candidate representation. The infimum need not be attained if the feasible set is not
closed (for example, strict tolerances or an open action domain). The algorithm must
report a policy at a declared budget and margin, never an alleged minimizer. For every
`eta>0`, if a policy of cost `<V+eta` has robust margin `gamma`, refinement eventually
finds a certificate when its enclosure error is below `gamma`. Conversely, certified
disjointness for all cells covering every policy of cost at most `b` proves `V>=b`. Equality can still occur when the infimum is not attained; a strict lower bound needs an additional separated argument. Exact optimality is asserted only
when closedness/compactness proves attainment.

## Complexity (only for this explicit policy class)

Let `N` be the number of policy cells after refinement, `M` the number of persistent
uncertainty boxes, and `K` the fixed number of time intervals. A straightforward
enumeration costs `O(N M K)` validated propagation calls, multiplied by the certified
integrator cost. For a uniform `d`-dimensional action grid of width `rho`, `N` scales as
`O(rho^-d)` (with separate factors for switching times and reset words); the state
dimension enters the enclosure cost and any state partition. This is not a complexity
claim for arbitrary observation-feedback policies, adaptive experiments, or infinite
dimensional systems.

For bounded sentences over Type-2 computable real functions, Theorem 2.1 of Gao, Kong, Chen and Clarke returns either that the delta-weakened formula is true or that the original formula is false. Where both hold, either answer is permitted. The positive relaxed answer need not satisfy the original constraints. See the [primary theorem, page 3](https://arxiv.org/pdf/1404.7171). This bounded hybrid reachability framework
does not yield universal decidability for unbounded, discontinuous, or misspecified
physical systems.

## Where weaker assumptions fail

* Treating an inner reachable union as if it were a robust policy certificate confuses
  `exists p forall theta` with `forall theta exists p`; different models may require
  different policies.
* Allowing `theta` to vary independently at each time step enlarges the trajectory set relative to a persistent parameter. This can make bounds unnecessarily conservative; it does not omit persistent trajectories when all parameter values remain admitted. Claiming exact reachability or completeness for the original persistent model needs the temporal coupling.
* Without a certified modulus/enclosure, corner samples do not bound an action cell; a
  narrow excursion can invalidate both positive and negative conclusions.
* Without finite candidate coverage, “all cells failed” proves nothing about omitted
  policies, and an infimum may be approached by an unenumerated sequence.
* Local certificates do not compose if shared baths, memory, timing, or conserved
  quantities are absent from the declared interface.
* Approximate simulation lacking alternating action refinement transfers observations,
  not executable control, so it cannot certify a prepared physical policy.

The constructive target is consequently modest but useful: an executable bounded policy
with an independently checkable outer-trace inclusion, or a complete budget-covering
separating witness, with `unknown` whenever the declared margin and model class do not
resolve the question.

## Primary sources

* G. Pola, A. Girard, P. Tabuada, “Approximately bisimilar symbolic models for
  nonlinear control systems,” https://arxiv.org/abs/0706.0246
* G. Pola, P. Tabuada, “Symbolic Models for Nonlinear Control Systems: Alternating
  Approximate Bisimulations,” https://arxiv.org/abs/0707.4205
* S. Gao, S. Kong, W. Chen, E. Clarke, “Delta-Complete Analysis for Bounded
  Reachability of Hybrid Systems,” https://arxiv.org/abs/1404.7171
