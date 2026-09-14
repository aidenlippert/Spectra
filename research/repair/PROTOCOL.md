# Representation repair from existing thermal evidence

Mission and bridge are unchanged. This pass constructs a missing predictive state
from observable records. It does not add another conditional architecture.

## Fixed workload and evidence boundary

Use the existing thermal equations and terminal requirements in
`experiments/mission_examples.py`: cold >=1/10 and hot <=3/10. Preserve all old
files and results. The original saved thermal examples contain nine reset,
two-control, endpoint observations. Consume all nine as initial evidence.

The learner receives controls and the two declared output ports, with exact
intervals as recorded in this mathematical workload. The adapter obtains these
ports from the existing returned task quantities, not from the hidden trajectory.
It receives no state matrix, hidden history, named missing state or candidate
answer list. The prior is deterministic, zero-reset, no-feedthrough, discrete-time
LTI behavior with one input, two output ports and dimension at most 4. This is an
explicit model-class assumption, not an inferred fact about real hardware.

Initial representation: current port 0 only. Analyze both discarded information
in the existing history and nonidentifiability of longer responses. Existing
two-step data identify the first two Markov vectors but not the whole future law.
If insufficient, request a reset impulse followed by zeros through 8 steps,
recording both outputs. Count reset, control steps and all scalar observations.
This provides the standard 2D Markov-data bound, not a claim of optimal experiment
selection. One pulse followed by zeros is already admissible in the old workload.

The learning executable may consume only its data JSON and use exact linear
algebra. Construct Hankel rows ordered by output port then future lag, select
independent rows/columns by exact rank, and derive predictive update, input and
output matrices. No physical matrix or intended rank is supplied. This is a
conventional Ho–Kalman/predictive-state baseline, not a novel algorithm or a
learned improvement over ARX/CLUE.

## Fixed consequential decision

Permit preparation histories followed by two operating steps under the same
equations, controls in [0,1], same final temperature constraints, and cost equal
to the sum of controls. Extending the history changes the requested task horizon
explicitly; the original two-step benchmark remains unchanged.

Check preparations (1,1/4) and (0,3/4), which have equal current hot output.
Search two operating controls from {0,1/4,1/2,3/4,1}, sorted by total control
then lexicographically, for a feasible continuation after preparation (1,1/4).
Compare the old representation/model-compatibility refusal with the repaired
prediction. Retain both cases and all policy failures, not only the passing one.
No claim of optimal preparation over all alternatives follows.

Before running acquisition, the first draft's (1,0) preparation was replaced:
its zero-input continuation is a prefix of the proposed training impulse and
would expose the policy answer directly. The revised preparation has a second
nonzero control; no four-step realization of it is in the initial records or
the impulse experiment. Thresholds, equations and control bounds are unchanged.

## Validation frozen before constructing the learner

* Independent checker uses the original equations to derive a linear map from
  physical to learned coordinates and checks exact intertwining, input, output
  and reset identities. This proves the supplied-model result for every finite
  admissible input history, not just sampled traces.
* Check every word over {0,1/2,1} of length 0 through 6 against the existing
  evaluator (1,093 words). These are fixed regression observations, not a claim
  of an unknown physical distribution or a novel cross-family discovery.
* Check deterministic feedback using only the preserved ports and input bounds.
  The algebraic equality covers any causal, well-defined port-feedback rule
  that leaves the component dynamics unchanged. It does not authorize an unknown
  reaction–thermal coupling; document that missing interface separately.
* Initialization is zero after the declared reset, or an exact online update
  through a recorded preparation history. Unrecorded history is refused unless
  a separately accounted observable initialization procedure is possible.
* Exact numerical/representation error is zero only for this exact model.
  Physical-model uncertainty, measurement noise, initialization error and any
  new coupling remain separate. Any uncertain observation is refused by this
  exact learner; it is not rounded into an exact system.

Success discharges the need to supply a state description and its update rule
for this admitted input-output family. It does not discharge finite-order LTI
model applicability, sensor validity, physical calibration, arbitrary resets,
nonlinear chemistry or global synthesis.
