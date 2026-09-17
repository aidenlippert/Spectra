# Audit of the dimer response diagnostic

## Preconditions that must be stated

The retained-weight identity is correct for any normalized eigenvector
`H\psi=E_0\psi`, provided `u=P\psi`, `v=Q\psi`, and `0<p=\|u\|^2<1`:

\[
 \langle v,(QHQ-E_0)v\rangle
 =\langle u,(PHP-E_0)u\rangle.
\]

The Schur inverse additionally requires `E_0\notin\operatorname{spec}(QHQ)`;
for the usual ground-state response one should state the stronger premise
`E_0<\lambda_{\min}(QHQ)`.  A ground eigenvalue degeneracy can invalidate a
claim that a selected ground vector has the relevant overlap or that the
inverse exists. For the product-dimer family the ground state is unique and
has strictly positive P weight. Hence it cannot belong to Q, and in finite
dimension the restriction to Q has strictly higher minimum energy. Merely
knowing E0<0 and PHP=0 would not prove this for a general Hamiltonian.

The response orientation in the note is now consistent: `B=QHP` maps retained
vectors to eliminated vectors and
`\Sigma=B^*(C-E)^{-1}B`.  From the projected equation,
`(C-E_0)^{-1}B u=-v`, hence the operator-norm lower bound is valid whenever
`u\ne0` and the inverse exists.

## Exact dimer formulas

For the half-filled two-site Hubbard singlet block
`[[0,-2],[-2,8]]`, the ground energy is `4-2\sqrt5`.  Its no-doublon weight is

\[
 w=\frac12\left(1+\frac{2}{\sqrt5}\right)
   =\frac{5+2\sqrt5}{10}.
\]

For `L` tensor-product dimers, `E_0=L(4-2\sqrt5)` and `p=w^L` exactly.
The current `dimer_family.py` uses rational outward intervals for these
radicals.  Its undressed gap bound is consistent with
`(-E_0)p/(1-p)`, using conservative interval directions.

The dressed diagnostic needs separate interpretation: a rational rotation
changes `PHP`, so the bare formula `PHP=0` no longer applies. The code's
`local_a` is the exact expectation in the normalized projection of the local
ground state onto the dressed retained space. Its product expectation is L
times that number; it is not a bound on every retained spin state.
It is not evidence that the dressed projector has the same Hubbard
no-doublon coercivity.

## Coercivity hypotheses

The sufficient-condition note should always enforce both
`a_*\ge0` and `U-a_*\ge0`, in addition to `U-a_*-b_*-E_max>0`.  Without
`a_*\ge0`, the stated local inequality does not have the intended monotone
interpretation; without `U-a_*\ge0`, the step
`(U-a_*)D\succeq(U-a_*)Q` reverses or becomes invalid.  Also, the relative
Neumann estimate requires `A` positive on the declared Q-sector and
`\rho=\|A^{-1/2}RA^{-1/2}\|<1`; those are premises, not consequences of a
local hopping norm.

## Scope

The dimer family rules out uniform conditioning at the exact ground energy
for one globally defined bare no-doublon elimination. It does not obstruct
finite-tolerance energy calculation, local, dressed, or adaptive
projectors, nor does it establish a lower bound on approximation degree at a
fixed tolerance.  Any stronger claim would require a separate theorem.
