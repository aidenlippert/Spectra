> Proposal status: unexecuted, not an established V5 result. The initial full-matrix parameter-count advantage was rejected; equally informed baselines can exploit the admitted operator family. Later qualifications are essential.

# V5: black-box operator discovery with a bounded compounding test

Use a two-node dissipative linear network with hidden coordinates. An episode has
state (x_t\in\mathbb R^2), intervention (u_t\), and observed output

\[
x_{t+1}=A x_t+B u_t+\xi_t,
\qquad y_t=Cx_t+\nu_t,
\]

where the learner sees only interventions and noisy observations. The physical
coordinates are hidden by an unknown invertible mixing (R): the supplied
interface exposes (z=Rx). Candidate mechanisms are sparse in the physical
coordinates: one episode has a diagonal damping operator, another has a single
off-diagonal coupling. No coordinate labels, support labels, or selected model
are supplied.

## First discovery

Task A is to identify an intervention-stable one-dimensional invariant mode from
short impulse-response experiments. The learner proposes a linear observable
(w^Tz), then actively chooses impulses that maximize disagreement among
candidate projectors. A held-out intervention tests whether the proposed mode
remains invariant when the input channel changes. Retain the mode only with a
confidence interval on the residual dynamics below ε; otherwise abstain.

The retained object is a reusable operator relation (a projector and validity
conditions), not a name such as “mode 1.” The certificate includes the observed
intervention class, noise bound, and residual bound.

## Second discovery and compounding criterion

Task B introduces a new coupling in the same hidden coordinate system. It starts
with the certified projector from Task A and tests only the orthogonal residual
subspace before searching the full two-dimensional operator family. A reset
solver repeats the full identification experiment. Both use identical shot,
intervention, precision, and computation budgets.

The strict first-order criterion is

\[
C_B^{\mathrm{warm}}<C_B^{\mathrm{reset}},
\]

where (C) counts physical interventions and observations, not only optimizer
steps. This can be ordinary reuse. A stronger compounding criterion requires a
third task C whose prediction error or intervention budget satisfies

\[
G_{AB}>G_A+G_B,
\]

with (G_S) measured against a frozen baseline under matched total resources.
The inequality must be evaluated on new couplings and interventions, not the
training trajectories. A conventional stateful solver with the same memory and
access is a mandatory baseline; if it achieves the same result, report a tie.

## Finite executable theorem

Restrict (R) to a finite grid with minimum singular value κ, restrict each
operator coefficient to a finite grid, and assume independent sub-Gaussian noise
with parameter σ. For every pair of hypotheses, include an impulse experiment
whose predicted output distributions differ by at least Δ in one measured
coordinate. Hoeffding/Gaussian concentration gives

\[
n\geq O\!\left(\frac{\sigma^2}{\Delta^2}
\log\frac{|\Theta|}{\delta}\right)
\]

shots for finite-family identification with confidence (1-δ). A projector
learned in Task A reduces Task B's hypothesis family from Θ to the residual
subfamily Θ\(_\perp\); the intervention saving is certified only if the
corresponding pairwise separation remains Δ under the new input channel.

If two hypotheses agree for every allowed intervention, they are one operational
equivalence class and must not be distinguished. If Δ is below noise, the
learner must abstain. This prevents the apparent compounding result from being
an artifact of an unmeasurable coordinate system.

## Negative controls and limits

Use a task where the first projector is invalidated by a changed coupling. The
correct output is certificate failure and a reset, not forced reuse. Use a third
task whose dynamics include a genuinely new latent mode; the old representation
may improve initialization but cannot certify that mode. Also compare against a
solver that stores all raw trajectories and reruns the same finite-grid search.

The construction demonstrates bounded cumulative operator reuse only when all
three conditions hold: a retained intervention-stable relation, a measured cost
reduction on the second discovery, and a superadditive held-out gain beyond the
strong stateful baseline. It does not establish general physical law discovery
or an advantage over an unrestricted emulator.
