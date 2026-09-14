# Repair implementation audit

The end-to-end run is reproducible and its learner is evidence-only. `parse_evidence`
consumes controls and exact output intervals, requires confirmed zero reset, and
does not receive the thermal matrix, hidden state, intended rank, or missing-state
name. The first two Markov vectors are identified from the initial records; the
rank-completion experiment is requested only after the exact moment test reports
`needs_evidence`.

The Hankel construction is consistently port-major and lag-major:
`H[(port,lag), column] = g[column+lag][port]`. The selected independent rows and
columns produce a rank-two realization, and the independent `equation_check`
verifies output preservation, transition intertwining, input intertwining, and
reset. The 1,093 fixed words and the port-feedback loop are post-fit checks; they
do not feed back into learning or policy selection.

The revised preparations `(1,1/4)` and `(0,3/4)` are absent from both the initial
records and the eight-step impulse acquisition. The operating search uses only
the learned model, then evaluates the selected policy with the pre-existing
thermal equation. The second evidence-compatible model agrees on all initial
records but makes the alternative preparation infeasible, so the policy-status
change is consequential rather than a cosmetic latent-state change. The selected
continuation is feasible under the declared cold and hot thresholds.

The method is correctly labelled as a conventional exact finite-order LTI
predictive realization, not a novel algorithm: it is a small Ho–Kalman/predictive
state construction with exact RREF. The run also separates physical-model,
measurement, initialization, numerical, and coupling uncertainty and makes no
global-optimality or physical-validation claim.

One scope caveat remains: the alternative witness is deliberately an
evidence-compatible mathematical model, not a claim about a second physical
thermal apparatus. The result therefore establishes a data-identification and
policy-repair construction conditional on the declared LTI class and the
supplied thermal equations. It does not establish applicability to real hardware
or to an unapproved reaction–thermal coupling.

## Boundary audit

`verify_against_records` now closes the main self-consistency gap: it requires
full-rank supplied design data, exact reproduction of every observable record,
and agreement between the claimed rank and reconstructed Hankel rank. The
adversarial test rejects an internally consistent fabricated model when checked
against raw evidence.

The bridge handoff now uses the learned prediction for its point `Sample`, while
the pre-existing thermal equation is checked independently for equality. This
resolves the prior attribution issue: the bridge receipt is downstream of the
repaired learner, with the equation evaluation retained as a separate post-check.

Cost fields charge records, scalar observations, prefix evaluations, rank tests,
fixed-word checks, and feedback steps. The checking timer now starts before model
unpacking, verification, and evidence binding, so those checks are included in
the reported checking bucket. These remain mathematical/computational costs,
not laboratory costs.

Resource boundaries are explicit: `D<=8`, at most four ports, 128 records,
fit histories of length at most `2D`, replay histories of length 256, and bounded
exact rationals. Malformed evidence and record mappings now fail through explicit
`ValueError`/`TypeError` refusal paths, including the CLI boundary. All scientific
claims remain conditional on the finite-order LTI prior and exact reset/port
evidence.

The final report labels the hot-only representation as a diagnostic ablation and
the compatible countermodel as manually constructed independent analysis. Neither
is attributed to learner invention.
