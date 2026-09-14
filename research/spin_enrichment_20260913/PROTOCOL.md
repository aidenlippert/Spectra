# Continue the small coupled spin-completion subspace

Keep the frozen rational H6 Hamiltonian, ten spatial patterns, collective
tail, quadratic baseline, body-one number ideal, and coefficient-L1 residual
rule. Continue from the accepted 64-direction result in
results/spin_completion_20260913/campaign/adaptive/round_4_batch_2.
No preceding source or artifact is modified.

The continuation has a 360-second outer watchdog, including imports, exact
source and upper replay, model reconstruction, pricing, every candidate
solve, export, exact acceptance and writes. At most twelve enrichment rounds
and 256 selected combinations are allowed. Each round proposes up to two
new numerically independent negative directions per exact symmetry block.
Directions are rounded before solving. Test bundles adding one and two per
block, without exceeding the dimension of that block or the total cap.
The unchanged solver has a 40-second numerical cap. Do not start a candidate
with less than 60 seconds left; the outer watchdog is authoritative.

Continue with the candidate giving at least 10^-6 Ha exact lower gain and
the best gain divided by trial wall time plus the round's pricing cost.
All candidates, including unselected and failed ones, remain charged.
Keep the strongest accepted lower independently of the continuation path.
If an accepted interval reaches 1.6 mHa, stop immediately; do not test a
larger bundle merely to improve an already sufficient bound. This gives a
first sufficient sampled subspace, not a global minimum-size claim.

Replay the inherited certificate before using its lower for scoring. Bind
its saved directions and frame references to the exact certificate; treat
its numerical dual only as a proposal. Reuse the previous cached CAR maps,
pricing, coordinate transforms, exporter, SOS checker and tail checker.
Do not reuse a wavefunction or prior SOS factors as discovery directions.

After continuation, the selected proof's directions receive one separate-
contribution ablation under a 90-second outer watchdog and 40-second solver
cap. It inherits all direction discovery cost. No full-family solve is run.
Fresh standard-library replay covers the source, all accepted new intervals
and the ablation, against the same rational upper. Report every cost and
failure, proof size, direction coefficients, largest blocks and target result.

Prior causal descriptor discovery cost is 222.510846501 seconds, including
the original counterexample's full-frame source solve and its exact repair,
spin-completion diagnostic and initial adaptive search. Earlier full-family
timeouts and ablation remain in the preceding whole-pass ledger; they are
not erased or charged as free successes. Incremental continuation costs and
cumulative causal discovery costs must be distinguished.

If pricing exhausts negative directions before the accuracy target, this is
only a numerical diagnostic. A full-family insufficiency claim still needs
an exact feasible dual and a separately recorded bounded proof phase.
