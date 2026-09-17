# Local PSD-cone obstruction for full-rank local marginals

This is an exact structural limitation on certificates whose positive terms
have bounded spatial support.  It applies to a fixed ground state and does not
assert that every lower-bound family must be a sum of local PSD terms.

## Exact-equality theorem

Let `H` act on a finite lattice and let `|\psi>` be a normalized ground state
with energy `E_0`.  Suppose an exact local certificate has the form

\[
 H-E_0=\sum_{C\in\mathcal C} h_C,
 \qquad h_C\succeq0,
\]

where each `h_C` acts only on a cluster `C`.  Let
`\rho_C=\operatorname{Tr}_{\bar C}|\psi\rangle\langle\psi|`.  Taking the
ground-state expectation gives

\[
 0=\sum_C\operatorname{Tr}(\rho_C h_C).
\]

Every summand is nonnegative.  If `\rho_C` is positive definite on the local
Hilbert space on which `h_C` acts, then

\[
 \operatorname{Tr}(\rho_C h_C)
 \ge \lambda_{\min}(\rho_C)\operatorname{Tr}(h_C),
\]

so every `h_C=0`.  Hence a nonzero exact frustration-free decomposition of
`H-E_0` into bounded-support PSD terms requires every nonzero term to have
support in a cluster whose ground-state marginal has a kernel.  Generic
entangled states have full-rank small marginals, so exact local PSD equality
is then impossible except trivially.

This is a cone statement, independent of commutativity.  It does not apply to
terms with indefinite sign, nonlocal response terms, or a certificate that
uses identities and cancellations outside the local PSD cone.

## Quantitative approximate version

Suppose instead a proposed lower certificate is

\[
 H-L=\sum_C h_C+R,\qquad h_C\succeq0,
 \qquad \|R\|\le\varepsilon.
\]

For the true ground state,

\[
 E_0-L=\sum_C\operatorname{Tr}(\rho_C h_C)+\langle\psi|R|\psi\rangle,
\]

and therefore

\[
 E_0-L\ge
 \sum_C q_C\operatorname{Tr}(h_C)-\varepsilon,
 \qquad q_C=\lambda_{\min}(\rho_C).
\]

Thus a family of patches with certified `q_C\ge q>0` cannot approach exact
energy equality unless the total positive mass
`\sum_C\operatorname{Tr}(h_C)` tends to zero or the residual allowance is
large.  This is a genuine quantitative obstruction once `q_C` and a lower
bound on the positive mass are independently certified.

## How to obtain a coefficient witness

The last lower bound is not implied by a large Hamiltonian coefficient norm:
different local PSD terms can cancel in their sum after embedding.  To make it
rigorous, define the linear embedding `\Phi(h)=\sum_C\iota_C(h_C)` and restrict
to the declared PSD cone.  A valid witness is a linear functional `\ell` on
Hamiltonian coefficients satisfying

\[
 \ell(\Phi(h))\le K\sum_C\operatorname{Tr}(h_C)
 \quad\text{for all }h_C\succeq0,
\]

and an independently computed inequality
`\ell(H-L)\ge b>0`.  Then

\[
 \sum_C\operatorname{Tr}(h_C)\ge b/K,
\qquad
 E_0-L\ge q b/K-\varepsilon.
\]

The inequalities defining `\ell` can be checked locally as PSD inequalities,
so the witness need not enumerate the global many-body basis.  Symmetry and
fermionic identities must be quotiented before applying this test; otherwise
`\ell` may falsely detect a coefficient that is algebraically null.

## Relevance to an eight-site Hubbard test

For an H8 model, one can enumerate only the local operator algebra on patches
of at most four sites, compute exact rational lower bounds on the eigenvalues
of the relevant local marginals, and search for a coefficient witness for the
embedded patch map.  A positive witness would certify that *that specific
bounded-support PSD cone* cannot reach the desired width without a residual
or nonlocal/cross-patch terms.  A failed witness would establish nothing about
the cone's sufficiency or the optimizer.

This test should use the actual declared H8 ground-state or certified-state
marginals and charge/spin sector, with all trace and normalization conventions
explicit.  It must not call a numerical full-state result an exact ground
state.  The result would be a precise obstruction to a local PSD certificate
family, while leaving open collective Schur responses and adaptive patch
systems—the mechanisms a larger theorem would need.
