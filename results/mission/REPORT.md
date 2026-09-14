# Constructive prediction, synthesis and invention: first mission audit

The [standing mission](../../research/MISSION.md) now explicitly preserves general
mastery of matter as the objective. Compounding and faster solvers are supporting
methods. The V12 cost test is closed; it failed its declared confidence screen.

The present result is a **conditional constructive connection**, developed and
checked against counterexamples. It is not a novel universal theorem or a newly
demonstrated physical capability. The full derivation is
[From physical evidence to reachable, certified design](../../research/mission/CONSTRUCTIVE_BRIDGE.md).

## The mathematical connection

Optimize complete preparation-and-operation policies, not only final-state
descriptions. Uniform evidence bounds over a covering of the admissible policy
space give two independently necessary objects:

- a feasible executable candidate with objective upper bound U;
- a lower bound L over every not-excluded policy cell.

Then the candidate's regret over that policy class is at most U-L. Every cell
must be excluded before claiming infeasibility. The result tolerates adaptive
selection only when the bounds are valid simultaneously across selection and
stopping. It keeps measurement uncertainty, execution risk, numerical error and
physical-model discrepancy distinct.

The finite experimental construction exposes costs usually omitted from the
prediction-to-design argument. Its certificate gap is bounded by

\[
 \eta+4w_0+2\omega_0(r)+\kappa(s),\qquad
 s_j=2w_j+\omega_j(r).
\]

Here eta measures a robust feasible competitor's suboptimality, w the site
uncertainty, omega policy-to-response continuity, and kappa the value gained by
relaxing feasibility. A useful budget therefore needs both reachable feasibility
margin and quantitative constraint conditioning. Lipschitz smoothness alone
does not supply them. A uniform grid remains exponential in policy dimension;
this is the explicit conventional baseline to improve, not a hidden efficient
oracle.

For coupled systems the error budget needs component state/input sensitivities
and interface gains. Stable state prediction does not automatically bound failure
probabilities at a threshold: an additional probability-of-near-boundary term is
required. These obligations directly connect representation, preparation history,
coupling and physical reliability.

## Executable verification and its exact scope

[mission_policy_cover.py](../../experiments/mission_policy_cover.py) constructs
and rechecks rational policy-space covers and response bounds. It rejects cover
holes, contradictory data/moduli, and unsupported guarantees. Insufficient data
or an unresolved noise floor remains unresolved. The input intervals and physical
policy decoder remain external premises; the code does not certify their origin.

Two distinct supplied process equations check the constructive argument:

| Mathematical family | Constructed two-step control policy | Certified objective interval | Known optimum | Exact oracle calls |
|---|---|---|---:|---:|
| Conserved two-state reaction | (1/32, 15/16) | [25/32, 31/32] | 9/10 | 32 |
| Two-node thermal transport | (13/16, 1/8) | [3/4, 15/16] | 4/5 | 9 |

Both retain their preparation trajectories and satisfy the supplied equations'
terminal constraints. The reaction evolves
`x_next=x+u(1-x)/2-x/4`, starts at x=0, and requires final x>=9/20. Its final
conversion is at most half total control, proving cost>=9/10; policy (0,9/10)
attains that bound. The thermal model evolves
`hot_next=hot/2+cold/4+u/2`, `cold_next=hot/4+cold/2`, starting at zero. It requires
final cold>=1/10 and hot<=3/10. Final cold=u0/8 implies cost>=4/5; (4/5,0)
attains it. These independent analytic optima check the interval direction.

The observations are exact evaluations of stated equations, **not laboratory
measurements**. No rates, mechanisms, representation or physical applicability
were autonomously discovered. These examples check mathematics without turning
their tiny domain into the mission.

[Saved policies, trajectories and evidence](mathematical_examples.json) ·
[replay and preservation receipt](verification_receipt.json) ·
[test log](tests.log).

## What remains open

The central unresolved problem is constructing affordable physical enclosures
and reachable policy descriptions across independently specified families when
the right mechanisms, state variables and preparation histories are not supplied.
This turn neither removes that obstruction nor demonstrates cross-domain
physical invention. It makes the obligation and a baseline construction precise.

The next central attack must replace one costly premise with evidence: discover
an intervention-sufficient representation, prove/check its coupling and history
domain, and use it to reduce physical prediction or synthesis cost on another
family. A learned confidence radius, a hand-supplied mechanism, or another small
solver timing improvement cannot count as that result. Existing primary results
and the boundaries of their relevance are recorded in the
[source map](../../research/mission/primary_source_map.md).
