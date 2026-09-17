# Audit of the PSD-seed Chebyshev lower bound

The proposed implication is sound with two bookkeeping corrections.

Let `W\succeq0` be a normalized Hilbert-Schmidt ground eigenmatrix,
`\|W\|_2=1`, and let the amplitude dimension be `d`.  For
`v_0=I/\sqrt d`,

\[
 \langle W,v_0\rangle=\operatorname{Tr}(W)/\sqrt d\ge1/\sqrt d,
\]

because `Tr(W)\ge\|W\|_2=1`.  This uses the proven PSD-ground-eigenmatrix
property of the real attractive Hubbard amplitude map.

For

\[
A=\frac{b+\ell-2H}{b-\ell},
\]

the ground eigenvalue satisfies `E_0<L` exactly when its Chebyshev coordinate
`x_0=(b+\ell-2E_0)/(b-\ell)` exceeds

\[
x_*=(z+z^{-1})/2,
\qquad
L=(b+\ell)/2-(b-\ell)(z+z^{-1})/4.
\]

For `0<z<1`, `T_k(x_0)\ge z^{-k}/2`.  Unrolling the inexact recurrence gives
the ground component of `v_k` as `T_k(x_0)\alpha` plus error terms.  The
standard Chebyshev estimate

\[
\frac{|U_{k-i}(x_0)|}{T_k(x_0)}
\le\frac{2z^i}{1-z^2}
\]

yields the sufficient exclusion condition

\[
2z^k\|v_k\|_2+
\frac{2}{1-z^2}\sum_{i=1}^kz^i\eta_i
<\alpha,
\qquad \alpha\ge1/\sqrt d.
\]

If the initial seed is approximate, with `v_0=I/\sqrt d+r_0`, subtract
`\|r_0\|_2` from the right side (or include it as an `i=0` error with the
corresponding coefficient).  Strict inequality is required; equality does
not exclude the hypothesis.

No spectral upper bound `b\ge E_{\max}` is needed for this one-sided logical
implication: only the ground eigenvalue is tested, and `b>\ell` defines the
polynomial.  A valid declared `b` is still useful for interpreting the full
spectral interval and for stable recurrence construction.  The energy claim
is only `E_0\ge L` for the same finite amplitude Hamiltonian and sector whose
ground PSD property was proved.

This route can avoid constructing dense `C` and `D`, because it needs matrix
vector actions and certified vector-norm/recurrence residual bounds rather
than a global PSD residual.  It does not automatically make those actions
cheap: exact rational squared norms can grow rapidly, and every omitted
arithmetic or approximation error must enter the `\eta_i` ledger.  It also
relies on the PSD-ground theorem and the seed overlap, so it is not a generic
many-body lower bound or a proof of dimension-independent complexity.
