# Limits of fixed-depth dressing on repeated dimers

This note sharpens the product-dimer obstruction.  It concerns the conditioning
of one *global* retained projector, not the cost of every possible algorithm.

## Exact product statement

Let one dimer have normalized ground state `|g>` and let `P_1` be a local
retained projector with

\[
 w=\langle g|P_1|g\rangle\in(0,1).
\]

For `L` disconnected identical dimers, take the product ground state
`|G_L\rangle=|g\rangle^{\otimes L}` and the product retained projector
`P_L=P_1^{\otimes L}`.  Then exactly

\[
 \langle G_L|P_L|G_L\rangle=w^L.
\]

More generally, for a product of local dressed projectors `P'_j` with local
weights `w_j\le1`,

\[
 \|P'_1\otimes\cdots\otimes P'_L|G_L\rangle\|^2
 =\prod_{j=1}^L w_j.
\]

If a fixed local construction has `w_j\le1-\epsilon` on a positive density
of dimers, the global retained weight is at most
`(1-\epsilon)^{cL}`, hence exponentially small.  A fixed-depth circuit that
acts independently inside each disconnected dimer can avoid this only by
making the local retained weight tend to one; a local unitary that maps the
local ground state exactly into the retained subspace does so with weight one.
There is no contradiction: that is an exact local dressing, not a global
complexity lower bound.

## What the statement does not establish

The argument does not prove that every fixed-depth circuit has exponentially
small overlap.  A depth-one circuit containing the exact local dimer rotations
already gives weight one.  Nor does small overlap imply large cost at fixed
absolute energy tolerance.  The relative residual metric
`M=I+X^*X` can absorb large response amplitudes, and a product system can be
certified by composing local responses without ever forming the global
projector.

For weakly coupled dimers, a product formula is no longer exact.  If the
coupling is treated perturbatively, the relevant control quantity is a
certified inter-dimer response norm divided by the local excitation gap.  A
fixed-depth dressing with local error `\epsilon` can also accumulate a global
state-norm error of order `L\epsilon`; that accumulation is a limitation of a
global state-norm guarantee, not necessarily of an energy-density or local
observable guarantee.

## A route that bypasses the global-overlap problem

A viable theorem should use local additive certificates.  For disconnected
regions, prove local response bounds

\[
 \Sigma_j^-(E)\preceq\Sigma_j(E)\preceq\Sigma_j^+(E),
\]

and compose energies or energy densities with an additive error ledger.  For
coupled regions, retain only boundary channels and certify the inter-region
response in the relative metric.  The required measurable condition is then a
uniform bound on each boundary channel/rank and on the sum of coupling-error
terms per region, rather than a lower bound on one global retained overlap.

This distinction is essential: the dimer family disproves the sufficiency of a
global bare-sector weight or a fixed local charge gap, but it leaves open a
local, dressed, composable response algorithm.  Any claim of impossibility
would need a lower bound on that compositional boundary complexity, not merely
on `\|P_LG_L\|`.
