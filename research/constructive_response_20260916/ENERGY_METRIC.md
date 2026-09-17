# Response error as an energy error

Let

\[
 S(E)=A_P-EI-B^*(C-EI)^{-1}B,
 \qquad C-EI\succ0
\]

on a closed window (I=[e_-,e_+]). Writing (X(E)=(C-EI)^{-1}B), direct
differentiation gives the exact operator identity

\[
 S'(E)=-M(E),\qquad M(E)=I+X(E)^*X(E)\succeq I.
\]

This is useful because it converts a response error into a displacement of an
energy crossing without assuming an exact eigenvector.

## Bracketing theorem

Assume uniform certified bounds

\[
 mI\preceq M(E)\preceq MI \quad (E\in I),
\]

and a constructible Hermitian approximation \(\widetilde S(E)\) satisfying

\[
\|S(E)-\widetilde S(E)\|\le\delta \quad(E\in I).
\]

If \(\widetilde S(\widehat E)\) has a zero lowest eigenvalue and

\([\widehat E-\delta/m,\widehat E+\delta/m]\subset I\),

then the exact crossing (E_0), whenever it exists in this interval, obeys

\[
\boxed{|E_0-\widehat E|\le\delta/m.}
\]

The proof is order-theoretic. At \(\widehat E\), positivity of the
approximation gives \(S(\widehat E)\succeq-\delta I\). A normalized vector
in the kernel of the approximation has Rayleigh quotient at most \(\delta\)
for the exact Schur operator. There is no general operator upper bound
\(S(\widehat E)\preceq\delta I\). Integrating
\(S'=-M) gives, for (h\ge0),

\[
 S(\widehat E-h)\succeq S(\widehat E)+mhI,
 \qquad
 S(\widehat E+h)\preceq S(\widehat E)-mhI.
\]

Therefore the first operator is PSD and the second has a nonpositive lowest
Rayleigh quotient at (h=\delta/m). Since (C-E\succ0), Schur complement
inertia identifies the PSD side with an energy lower bound and the negative
side with an energy upper bound. This argument uses no common eigenbasis and
no exact eigenstate. Strict signs or an endpoint convention handle equality.

More generally, if an interval response certificate gives
\(L(E)\preceq B^*(C-E)^{-1}B\preceq U(E)\), set
\(S_- = A_P-E-U\) and \(S_+=A_P-E-L\). Checking (S_-(e_- )\succeq0)
and a nonzero vector with \(v^*S_+(e_+)v\le0\) supplies a rigorous bracket
directly. Positivity of \(S_+\) alone is not an energy upper certificate. The metric result
above quantifies how much the bracket can move when a trial crossing is used.

## Constructible metric from a trial response

If \(\widehat X\) is available and (R=B-(C-E)\widehat X), then
\(X=\widehat X+Z), (Z=(C-E)^{-1}R). A certified gap
\(C-E\succeq\alpha I) gives \(Z^*Z\preceq\alpha^{-2}R^*R\). Hence

\[
 I+\tfrac12\widehat X^*\widehat X-\alpha^{-2}R^*R
 \preceq M(E)
 \preceq I+2\widehat X^*\widehat X+2\alpha^{-2}R^*R.
\]

These are directly checkable local/operator inequalities if the trial action,
residual contraction, and gap certificate are exact. A relative preconditioner
can replace the scalar gap and sharpen the (Z^*Z) term. This construction
does not require forming (X), (C^{-1}), or a determinant matrix.

## What this does and does not establish

The theorem supplies a clean observable target: certify a uniform response
error \(\delta\), then report the induced energy allowance \(\delta/m\).
It prevents a large response norm or small wavefunction overlap from being
mistaken for an energy failure. It does not prove that \(\delta\) stays small
under repeated many-body elimination. The unproved part remains the
constructive, non-enumerating representation of the residual and its metric,
with controlled growth through interacting levels.
