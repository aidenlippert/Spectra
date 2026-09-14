# Review: finite-horizon policy-cover stopping theorem

The proposed rule is sound only as a conditional statistical decision procedure. It
does not by itself certify a physical design, model validity, or global optimality.
Let (E) be one simultaneous confidence event on which every reported enclosure
used by the run is valid for every policy/cell queried (including adaptively chosen
queries). All conclusions below must be explicitly conditioned on (E), with

\[
\Pr(E)\ge 1-\delta.
\]

## Quantifiers that must be explicit

The policy class Π must be finite, or have a finite, computable cover with a stated
metric and a proved continuity/modulus relation from cover cells to objective and
constraints. A cover generated after seeing data is valid only if the confidence
construction is uniform over the entire possible adaptive transcript; ordinary
pointwise intervals do not suffice. If the number of refinements is unbounded,
confidence spending (for example, a summable α schedule) or an anytime-valid
construction is required. Reusing observations in later intervals must preserve
the advertised simultaneous coverage.

For each policy (J) and constraint (g), say whether the enclosure is for the
true physical quantity, a conditional expectation, or a finite-horizon execution
probability. The latter needs a declared randomness model and a reachability/event
definition. A simulator interval is not a physical guarantee unless discrepancy
between simulator and experiment is bounded on the queried domain. Adaptive
experiments create selection bias unless the interval method accounts for the
selection rule.

The objective and constraints need a common direction. With minimization, a
feasible candidate is certified by (U_g\le0), while an outer cell can be excluded
only when its *best possible* constraint is positive (a valid lower bound), or when
its objective lower bound already exceeds the current feasible upper bound. A cell
whose constraint lower bound is nonpositive remains possible even if its point
estimate looks infeasible. For maximization all inequalities reverse.

## Regret-gap statement

On (E), let (U_{best}) be the objective upper bound of an actually feasible
candidate and (L_{best}) the lower bound over every not-excluded cell, including
the candidate's own cell where appropriate. Then the true optimum over the covered
domain lies between these values and the candidate regret is at most

\[
U_{best}-L_{best}.
\]

This requires the cell lower bounds to be valid for the infimum over the whole cell,
not merely for its representative. If a feasible candidate is absent,
`U_best` is undefined and the run must report “no certified feasible policy,” not a
regret result. If all cells are excluded, the conclusion is only “no policy in the
specified covered class satisfies the declared bounds,” subject to (E); it is not
an impossibility result for physical policies outside the class, outside the cover,
or under a misspecified model.

## Finite termination

Finite stopping follows from a finite initial cover, finite work per cell, and a
rule that either terminates, refines a cell, or excludes it. To guarantee eventual
termination with a regret target ε, one additionally needs a robust near-optimal
competitor: a policy/cell with true feasibility margin at least ρ>0 and objective
within the target, plus a known refinement modulus so that sufficiently small cells
have uncertainty widths below ρ and the remaining objective gap. Without a strict
feasibility margin, a boundary-feasible optimum may require infinite refinement.
Likewise, if confidence widths have no convergence guarantee, refinement need not
decide any cell. A finite policy class alone does not fix either problem when the
reported enclosures remain statistically noisy.

The robust competitor assumption must quantify reachability: the policy must be
executable under the declared finite horizon, resource limits, observation model,
and admissible interventions. If reachability is only hypothesized, the theorem
certifies a conditional search result, not construction of the physical system.

## Model misspecification and probability bookkeeping

Separate aleatoric execution randomness, measurement noise, numerical error, and
model discrepancy. Confidence intervals for the first two do not cover the last.
A valid theorem needs either a physical calibration bound, a discrepancy interval,
or an explicit assumption that the model is correct on the covered domain. If the
model assumption fails, the result can retain mathematical validity only for the
model, not for matter. The failure probability must include every adaptive stage,
cover-selection event, and any external validation decision; reporting a nominal
per-query confidence is insufficient.

The finite-horizon qualifier also matters: a guarantee for trajectories through
time (T) does not imply safety, stability, or property retention after (T), nor
successful synthesis beyond the observed preparation path.

## What the theorem can and cannot claim

With the conditions above, the rule gives a useful certificate: with probability at
least (1-δ), a returned candidate is feasible and has bounded regret over the
declared covered policy class, or every covered cell is ruled out. It does not prove
that the best physical policy was found, that the cover contains all reachable
materials, that an optimizer's numerical answer is globally optimal, or that a
model-discrepancy-free interval transfers to an experiment. Those stronger claims
need a separate reachability theorem and physical validation protocol.

## Audit of `CONSTRUCTIVE_BRIDGE.md`

The deterministic certificate is correct under its strong simultaneous cell-wise
enclosure quantifier. The finite-cover calculation, however, has an endpoint
counting error. If a net site is within (r) of (d_\eta), a site interval has
halfwidth (w_j), and cell expansion costs \(\omega_j(r)\), then

\[
U_j(C_i)\le g_j(d_\eta)+w_j+2\omega_j(r).
\]

Consequently the robust margin is not the displayed
\(s_j=2w_j+\omega_j(r)\); under the stated site-interval-plus-cell-expansion
construction it must be at least (w_j+2\omega_j(r)). The draft needs to define
whether (L_j(C),U_j(C)) are already cell-wide endpoints. Mixing those endpoints
with site intervals changes the count again.

The same ambiguity invalidates the displayed
\(\eta+4w_0+2\omega_0(r)+\kappa(s)\) bound as written. With direct cell endpoints,
the candidate upper error is (w_0+2\omega_0(r)), while a surviving-cell lower
endpoint is (J(d_k)-w_0-\omega_0(r)) before any additional representative-to-policy
comparison. If that comparison is required on both sides, another modulus term is
incurred. A corrected proof must fix one endpoint convention and derive (s) and
the objective gap from it; the current proof does not establish the stated formula.

The finite grid also requires a known bias floor, noise scale, finite execution cost,
and a rule for failed preparations. “All cells excluded” remains relative to the
decoded cover and simultaneous-event assumptions, never an impossibility result for
uncovered reachable chemistry or a misspecified physical model.

### Correction after notation clarification

The preceding endpoint objection does not apply if the theorem's intended objects
are distinct: a point candidate upper bound at the sampled center, and a cell lower
bound obtained by expanding the site's lower endpoint. In that convention, an
interval [μ−w,μ+w] containing (f) gives a point upper error (2w). For a
candidate center within (r) of (d_\eta), this yields
\(g(d_i)+2w\le g(d_\eta)+\omega(r)+2w\), so the robust condition
\(\rho\ge2w+\omega(r)) is correct. Likewise a surviving cell's lower endpoint
\(\mu_k-w-\omega(r)\le0\) implies the center value is at most that endpoint
plus (2w+\omega(r)), giving the stated relaxed-set parameter. The corresponding
candidate objective and cell objective errors sum to
\(4w_0+2\omega_0(r)), so the displayed gap follows under this explicit
point-versus-cell convention. The audit should therefore require that convention
to be stated clearly, rather than reject the formula.
