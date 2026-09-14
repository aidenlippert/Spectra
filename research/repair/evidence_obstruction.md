# Evidence obstruction: isolated reaction/thermal records do not identify coupling

This note audits the existing mathematical examples in
`experiments/mission_examples.py`. The file supplies two separate deterministic
input-output maps, `reaction(policy)` and `thermal(policy)`, and exact records
for each map. It supplies no law for their joint execution, shared resource,
interference, or cross-state. Therefore a coupled-process decision is not
identified by those records, even when both isolated equations are known
exactly.

## Coupling extension theorem

Let (R:U_R\to Y_R) and (T:U_T\to Y_T) denote the supplied reaction and
thermal maps, including their histories. Let (D_R\subseteq Y_R) and
(D_T\subseteq Y_T) be the target sets already present in the workload. A
joint policy is (u=(u_R,u_T)). The isolated evidence constrains only the
marginals

\[
 y_R=R(u_R),\qquad y_T=T(u_T).
\]

It does not constrain a joint feasibility variable (c(u_R,u_T)), or any
cross-effect that is absent from both isolated interfaces.

For any selected existing target pair (A_R\subseteq U_R), (A_T\subseteq
U_T), define two logical coupling extensions:

\[
 \mathcal W_+(u_R,u_T):
   [R(u_R)\in D_R]\land[T(u_T)\in D_T],
\]

and

\[
 \mathcal W_-(u_R,u_T):
   [R(u_R)\in D_R]\land[T(u_T)\in D_T]\land
   [ (u_R,u_T)\notin A_R\times A_T].
\]

Choose (A_R,A_T) to contain the already proposed target policy (or any
existing target policies under comparison). Both worlds have exactly the same
reaction-only and thermal-only answers for every input and every history. They
therefore agree on every isolated observation, including all values generated
by `mission_examples.py`. Yet the existing target policy is feasible in
(mathcal W_+) and infeasible in (mathcal W_-). No new threshold,
benchmark, or tuned target is introduced: the only predicate used is the
existing target predicate and the already proposed policy.

This is deliberately a logical extension, not a claimed physical coupling law:
the supplied examples declare no admissible coupling class from which a
mechanistic extension could be selected. If a physical class is supplied later,
the same nonidentification question must be rerun inside that class.

### Information-theoretic conclusion

Let (E) be the complete isolated evidence. The two extensions satisfy

\[
 P_+(E)=P_-(E)
\]

but have opposite truth values for the existing joint-feasibility claim. Thus
that claim is not a function of (E). Any decision rule based only on (E)
returns the same answer in both worlds and is wrong in at least one. This is
an identification obstruction, not a claim that the isolated equations are
wrong. A coupling certificate requires an interface law, a shared-state prior,
or a discriminating joint experiment.

## Finite records hide delayed higher-order state

The same issue occurs even for one subsystem when the record has finite input
and output horizon. Fix any finite horizon (T), and consider discrete-time
linear systems with scalar input (u_t) and output (y_t). For each integer
(m>T), define the delay-line system

\[
 x_{t+1}=S_mx_t+e_mu_t,\qquad y_t=e_1^\top x_t,
\]

where (x_t\in\mathbb R^m), (S_me_i=e_{i-1}) for (i>1), (S_me_1=0), and
the input enters through (B=e_m). Starting from (x_0=0), its input-output
relation is

\[
 y_t= u_{t-m} \quad (t\ge m),
 \qquad y_t=0\quad(0\le t<m).
\]

Compare it with the zero system (y_t\equiv0). For **every** input sequence,
the two systems produce identical outputs through time (T), because
(m>T). Hence no finite input-output record of horizon (T), even one
containing exact outputs under adaptively selected inputs, can distinguish the
zero system from all higher-order delayed systems. After time (m), the same
systems can respond differently to a preparation pulse, so a policy justified
from the finite record can fail when the delayed state becomes visible.

This is a discrete-time construction. It is sufficient for the information
obstruction; no exact finite-horizon delay claim for finite-dimensional analytic
continuous-time LTI systems is made here.

If the admissible state dimension is known to satisfy (m\le d_{max}), a
record long enough to expose every allowed delay can make this particular
obstruction finite. Without a maximum-dimension/order prior (or an equivalent
finite complexity bound, known realization class, or persistent experiment
that identifies it), there is no finite horizon that rules out all delayed
higher-order states. A fitted low-order model is therefore not evidence that
the omitted memory is absent.

## What a constructive repair must add

The existing examples can support a coupled policy only after one of the
following is supplied and checked:

1. an explicit coupling/interface law with a stated admissible class and
   bounded model discrepancy;
2. a prior that bounds shared-state dimension/order (plus an identifying input
   experiment and state initialization rule); or
3. a joint experiment that distinguishes the competing coupling extensions and
   probes the relevant delayed state before the policy is certified.

These are alternative information sources, not assumptions that can be inferred
from the two isolated functions. The exact isolated reaction and thermal
equations remain useful baselines, but they cannot by themselves justify the
joint feasibility claim requested by the mission.
