# A posteriori compact response bounds

This note isolates the strongest result currently justified by the response
idea. It is a theorem about a supplied operator action and a certified metric;
it is not yet a complexity theorem and it does not prove that a compact trial
response exists for generic many-body Hamiltonians.

## The theorem

Let (A(E)=C-EI) act on an eliminated sector (Q), and let (B=QHP) be
the coupling to a retained sector (P). Assume, on an energy window

\[
 E\in I=[E_-,E_+],\qquad 0<\alpha I_Q\preceq A(E)\preceq\beta I_Q.
\]

For a constructible trial response (X(E):P\to Q), define the exact residual

\[
 R(E)=B-A(E)X(E).
\]

Suppose compact operators (M_-(E),M_+(E)) are independently certified to
satisfy

\[
 0\preceq M_-(E)\preceq A(E)^{-1}\preceq M_+(E).
\]

Then the Schur response

\[
 \Sigma(E)=B^*A(E)^{-1}B
\]

has the exact, a posteriori Loewner enclosure

\[
\boxed{
 L_X(E)\preceq \Sigma(E)\preceq U_X(E),
}\]

where

\[
 L_X=2\operatorname{Re}(X^*B)-X^*A X+R^*M_-R,
 \qquad
 U_X=2\operatorname{Re}(X^*B)-X^*A X+R^*M_+R.
\]

Consequently the enclosure width is itself certified:

\[
 0\preceq U_X-L_X=R^*(M_+-M_-)R.
\]

If only (\|R\|\le\rho) and (A\succeq\alpha I) are available, one may
take (M_-=0,M_+=\alpha^{-1}I), giving

\[
 L_X\preceq\Sigma\preceq L_X+\alpha^{-1}R^*R
 \preceq L_X+\rho^2\alpha^{-1}I.
\]

The identity behind the theorem is the error representation. Writing
(x=A^{-1}B=X+A^{-1}R),

\[
 \Sigma=2\operatorname{Re}(X^*B)-X^*AX+R^*A^{-1}R.
\]

Thus no estimate of a cross term is needed: the cross term is absorbed exactly
by the quadratic variational identity. The result applies to matrices,
operators, or exact CAR polynomial actions, provided each displayed Loewner
inequality is actually checked in the relevant representation.

## How to certify the metric without enumerating determinants

The theorem does not make (M_\pm) free. A useful certificate must prove

\[
 M_+-A^{-1}\succeq0 \quad\Longleftrightarrow\quad
 A^{1/2}M_+A^{1/2}\succeq I,
\]

and (A^{-1}-M_-\succeq0) analogously. A sufficient, directly checkable
condition is an operator sandwich

\[
 mI\preceq A\preceq MI,
\]

because (M_-=M^{-1}I), (M_+=m^{-1}I) are valid. This scalar choice is
usually too loose, but it needs only a gap and norm certificate. A stronger
compact metric can be checked by proving (A^{1/2}M_+A^{1/2}-I\succeq0)
through an SOS/Gram identity, a block-local PSD decomposition, or a certified
MPO/automaton action bound. The checker must verify the identity and PSD
remainders exactly (or with a separately proved rounding allowance).

This gives a concrete non-enumerative acceptance test: construct (X,M_\pm),
construct the residual action (R=B-AX), and check the small/local PSD rules
plus the exact residual contraction. Merely fitting (M_\pm) to sampled
eigenvectors, or compressing a matrix after determinant enumeration, is not
such a certificate.

## Uniformity over an energy window

For (E\le E_+), (A(E)=C-EI\) is monotone decreasing. A single certified
gap (C-E_+I\succeq\alpha I) and norm bound (C-E_-I\preceq\beta I)
provide uniform scalar metrics. For a polynomial trial (X(E)), uniform
residual bounds follow from an exact polynomial identity or interval arithmetic
on the coefficients. If (M_\pm(E)) vary with (E), the checker must verify
their metric inequalities for the entire interval, for example by a positive
coefficient Bernstein/SOS certificate after mapping (I) to ([0,1]).

The resulting enclosure is uniform only if the residual and metric checks are
uniform. Checking a grid of energies is evidence for a diagnostic, not a
proof over the interval.

## Two elimination levels

Partition the first retained space again as (P=P_2\oplus Q_2). Apply the
theorem to the first block with (A_1=C_1-E), producing

\[
 \Sigma_1\in[L_1,U_1],\qquad
 S_1=A_1^{\rm ret}-\Sigma_1.
\]

Loewner order reverses under subtraction, so the exact second-level operator
obeys

\[
 A_1^{\rm ret}-U_1\preceq S_1\preceq A_1^{\rm ret}-L_1.
\]

Use the appropriate certified upper/lower enclosure as the input operator for
the second response, and retain every induced term in its residual action.
If its metric bounds are (M_{2,-}\preceq A_2^{-1}\preceq M_{2,+}), the same
formula gives a second enclosure. Substitution is valid because congruence
preserves Loewner order; the total reduced certificate is the composition of
the two checked PSD remainders. Widths add under the corresponding congruence
and subtraction maps, with amplification bounded explicitly by the operator
norms of the retained couplings. One must report this amplification; “the
errors add” without it is false in general.

## What is proved, conjectured, and blocked

## A constructible relative metric

Avoiding (A^{1/2}), let (Dsucc0) be a compact preconditioner and certify
the relative sandwich directly:

\[
 (1-\eta)D\preceq A\preceq(1+\eta)D,\qquad0\le\eta<1.
\]

If (D^{-1}) is available as an exact local/MPO action, inversion order gives

\[
 (1+\eta)^{-1}D^{-1}\preceq A^{-1}preceq(1-\eta)^{-1}D^{-1}.
\]

Thus a relative residual can be checked by the compact quadratic certificate

\[
 R^*D^{-1}R\preceq\omega I,
\]

which yields response width at most (2\eta\omega/(1-\eta^2)). This is a
real acceptance rule only when the sandwich and quadratic inequality are
verified by exact operator identities, local PSD blocks, or a separately
certified norm bound. Sampling or a post-enumeration fit is insufficient.

## Exact recursive rule

For a positive block matrix (M=\begin{bmatrix}A&B\\B^*&C\end{bmatrix}),
define (\mathcal S_C(M)=A-B C^{-1}B^*). If

\[
 0\prec C_-\preceq C\preceq C_+,
\]

then inverse order and congruence give the checked enclosure

\[
 A-B C_-^{-1}B^*\preceq\mathcal S_C(M)
 \preceq A-B C_+^{-1}B^*.
\]

At level one, retain the response bounds (L_1\preceq B^*C^{-1}B\preceq U_1).
The exact next block therefore has the *reversed* bounds

\[
 A_1-U_1\preceq S_1\preceq A_1-L_1.
\]

At level two, apply Schur monotonicity to the two explicitly constructed
endpoint operators, provided their eliminated blocks are positive. This
monotonicity holds even when cross blocks differ, because the Schur quadratic
form is the minimum of the full quadratic form over the eliminated variable.
Every induced term must be included.

For a quantitative error bound use the whole path M(t)=M_-+t Delta, in the
orientation M=[[A,B*],[B,C]]. If C(t)>=delta I and ||B(t)||<=g throughout,
the derivative of its Schur complement is V(t)* Delta V(t), where
V(t)=[I; -C(t)^(-1)B(t)]. Therefore

\[
 \|\mathcal S(M_+)-\mathcal S(M_-)\|
 \le [1+(g/\delta)^2] \|\Delta\|.
\]

This is a norm bound, not an operator comparison of matrices acting on
different spaces. The factor explains why errors need not simply add.
`composition_exact.py` checks actual successive eliminations of a
five-dimensional rational matrix, endpoints with cross-block perturbations,
and a positive-definite example where a unit eliminated-diagonal error changes
the Schur complement by eighteen units.

Proved here is the variational residual identity, the two-sided response bound,
the scalar-metric fallback, window uniformity under uniform certificates, and
the order rule needed for recursive substitution. These are algebraic facts.

The research conjecture is that molecular Hamiltonians of practical interest
often admit (X) and (M_\pm) represented by bounded-support or low-bond MPO/SOS
objects whose certified width remains useful under repeated elimination. No
such generic bound is known here. The existing Spectra work demonstrates
compact operator actions and individual response components, but does not yet
certify the terminal retained operator or end-to-end cost advantage.

The exact barrier is therefore precise: a compact response is valuable only
when (i) the residual is constructible without full determinant enumeration,
(ii) a nontrivial metric sandwich is independently certified in the same
representation, and (iii) recursive induced terms stay representable with
controlled bond/support growth. Failure of any one leaves a valid diagnostic,
not a world-level breakthrough.
