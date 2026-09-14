# Workload pair: a decision-relevant missing thermal state

The existing `mission_examples.thermal` equation is

\[
(h',c')=(h/2+c/4+u/2,\;h/4+c/2),
\]

starting from `(h,c)=(0,0)`, with the declared constraints `c >= 1/10` and
`h <= 3/10`. The initial hot-only representation is inadequate even before any
new physics is introduced.

Two reachable length-two input histories have the same current hot value:

| history | current `(hot,cold)` |
|---|---|
| `(0, 1/4)` | `(1/8, 0)` |
| `(1/2, 0)` | `(1/8, 1/16)` |

Thus a representation retaining only current hot identifies the histories. Apply
the same declared continuation `(0, 1/2, 0)` to both. Exact iteration gives:

| history | final `(hot,cold)` | policy status |
|---|---|---|
| `(0,1/4)` | `(39/256, 45/512)` | infeasible (`45/512 < 1/10`) |
| `(1/2,0)` | `(169/1024, 13/128)` | feasible (`13/128 >= 1/10`, hot below `3/10`) |

The final hot values are both below the hot limit. The feasibility difference is
therefore caused by the omitted cold memory, not by a threshold chosen after the
fact. A predictor or policy decoder using only current hot must return the same
decision for both histories and necessarily misclassifies at least one.

The repair is the explicit two-coordinate state `(hot,cold)`, initialized by the
two observed input transitions and updated by the already declared thermal
equation. This is a known conventional realization of the supplied equations,
not a discovery of new physics or a claim of autonomous representation search.
The constructive cost is two input transitions to initialize the state plus one
update per subsequent control; exact rational replay is constant work per step.
For an actual workload, the missing obligations remain state observability,
physical applicability of the two-node model, and validation under declared
couplings. No physical measurements are present in this example.
