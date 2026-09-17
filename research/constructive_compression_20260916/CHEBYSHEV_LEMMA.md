Exploratory note. The reviewed derivation and implemented scope are in [DERIVATION.md](DERIVATION.md).

# Inexact Chebyshev PSD witness

Let (H) be self-adjoint with certified spectral upper bound (H\preceq bI),
and choose (\ell<b). Define

\[
 A={b+\ell-2H\over b-\ell}.
\]

If the ground energy (E_0<\ell-\delta), then the ground eigenvalue of (A)
satisfies

\[
 x_0>1+{2\delta\over b-\ell}=\cosh\theta.
\]

Let (W\succeq0) be a Hilbert-Schmidt-normalized PSD ground eigenmatrix supplied by the
reflection-positive cone argument, and let (v_0=I_N/\sqrt d), where (d) is
the dimension of the one-spin sector. If

\[
 |\langle W,v_0\rangle|\ge d^{-1/2},
\]

construct inexact Chebyshev vectors

\[
 v_1=Av_0+r_1,\qquad
 v_j=2Av_{j-1}-v_{j-2}+r_j.
\]

Assume (\|r_j\|\le\eta_j), with the norm and residual convention fixed by
the checker. Projection onto (W) gives

\[
 T_k(x_0)\langle W,v_0\rangle
 =\langle W,v_k\rangle-
 \sum_{i=1}^kU_{k-i}(x_0)\langle W,r_i\rangle.
\]

Therefore a contradiction to (E_0<\ell-\delta) is certified whenever, with
(x_*=1+2\delta/(b-\ell)=\cosh\theta),

\[
 {\|v_k\|\over T_k(x_*)}+
 \sum_{i=1}^k {U_{k-i}(x_*)\eta_i\over T_k(x_*)}
 < {1\over\sqrt d}.
\]

The monotonicity used here is elementary: for (x>1), both (T_k(x)) and
(U_n(x)) are positive, and (U_{k-i}(x)/T_k(x)) decreases as (x) grows
for fixed indices in the displayed range. Thus replacing (x_0) by (x_*) is
conservative.

For exact rational-friendly checking, set (z=e^{-\theta}\in(0,1)), so
(x_*=(z+z^{-1})/2). Since

\[
 {1\over T_k(x_*)}\le2z^k,
 \qquad
 {U_{k-i}(x_*)\over T_k(x_*)}
 \le {2z^i\over1-z^2},
\]

the simpler sufficient test is

\[
\boxed{
 2z^k\|v_k\|+
 {2\over1-z^2}\sum_{i=1}^kz^i\eta_i
 <{1\over\sqrt d}.
}
\]

The corresponding certified lower bound is

\[
 E_0\ge {b+\ell\over2}-{b-\ell\over4}(z+z^{-1}).
\]

This is a spectral witness, not an eigenstate approximation: it needs only
the PSD overlap premise and certified recurrence residuals. It does not remove
the need to construct and verify the MPS operator actions, and it may fail
when the overlap (d^{-1/2}) is too small or residual accumulation dominates.
The lemma is conditional and should not be presented as a new general
many-body algorithm until an actual compact recurrence passes on the target.
