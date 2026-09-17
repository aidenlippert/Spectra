# Local coercivity for doublon elimination in the Hubbard model

Status: candidate sufficient condition, not an achieved theorem.

For a finite graph, write

\[
H=UD+T,\qquad D=\sum_i n_{i\uparrow}n_{i\downarrow},
\]

with (T=-t\sum_{(ij),\sigma}(c^\dagger_{i\sigma}c_{j\sigma}+h.c.)).  Let (P)
be the global no-doublon projector and (Q=1-P).  Holes and charge transfer
remain allowed; fragment particle numbers are not fixed.  The exact response is

\[
\Sigma(E)=PHQ(QHQ-E)^{-1}QHP.
\]

The required premise is a certified lower bound on (QHQ-E), not the claim that
U itself is a global gap: (QTQ) can move doublons and lower the energy.

## A local sufficient condition

Suppose finite clusters (C), with weights (w_C>0), have exact inequalities

\[
T_C\succeq-a_CD_C-b_CI_C,
\qquad D_C=\sum_{i\in C}n_{i\uparrow}n_{i\downarrow},
\]

where boundary hopping is assigned exactly once.  Overlap accounting then gives

\[
T\succeq-a_*D-b_*I,
\quad a_* = \max_i\sum_{C\ni i}w_Ca_C,
\quad b_* = \sum_Cw_Cb_C.
\]

This implication requires the explicit premise (U-a_*\ge0), because only then
does ( (U-a_*)D\succeq(U-a_*)Q) follow from (D\succeq Q) on the eliminated
sector.  Under that premise, (QHQ\succeq(U-a_*-b_*)Q).  For an energy window
(E\le E_{\max}), a checkable eliminated-sector gap is

\[
\Delta_Q=U-a_*-b_*-E_{\max}>0.
\]

If (QHQ-E=A+R), (A\succeq\Delta_QQ), and
\(ho=\|A^{-1/2}RA^{-1/2}\|<1), then the directly constructible response
\(\Sigma_0=PHQA^{-1}QHP) obeys

\[
\|\Sigma-\Sigma_0\|
\le \|PHQA^{-1/2}\|^2\frac{\rho}{1-\rho}.
\]

The (m)-term Neumann response has remainder at most
\(|PHQA^{-1/2}\|^2\rho^{m+1}/(1-\rho)).  This is the compact-response
criterion worth testing: local certificates establish (A) and ρ, while only
boundary couplings are retained explicitly.

The useful form is relative:

\[
|\langle\psi,T\psi\rangle|\le a\langle\psi,D\psi\rangle+b\|\psi\|^2,
\qquad a<U.
\]

A norm estimate (\|T\|=O(|V|t)), or (U-z t|V|), is generally useless in the
thermodynamic limit.  Doping and charge transfer are handled by retaining the
global particle sector and never imposing a particle number per fragment.

## Frustration and the obstruction

Next-neighbor hopping changes the cluster constants and boundary response.  All
signed hopping terms must be included; absolute-value row sums may certify a
bound but can make (a_*\ge U).  Failure of this criterion means the selected
cluster radius or response order is inadequate, not that the model is gapped or
gapless.

The naive claim (QHQ\succeq UQ) is false whenever a normalized doublon-hole
state has (\langle QTQ\rangle<0).  Local doublon costs alone therefore do not
certify an eliminated gap; adjacent doublon-hole motion belongs to (QTQ) and
must be charged in the local inequality.

This is related to, but stronger than, a formal Schrieffer--Wolff/(t/U)
expansion.  Known strong-coupling theory gives the (4t^2/U) exchange term under
an assumed separation.  The proposed new result would be an explicit local
coercivity certificate and composable response remainder, uniform over declared
doping/frustration conditions.  Reproducing the exchange identity is not new.

## Decisive bounded experiment

On (2\times2) and (2\times3) clusters at (U/t=8), certify rational local
lower bounds for (T_C+a_CD_C+b_CI_C) in half-filled and one-hole sectors, for
nearest- and next-neighbor hopping.  Optimize overlap weights subject to
(a_*<U), then evaluate Δ_Q and ρ.  A positive stable pair (Δ_Q,ρ<1) is a
reusable construction; an exact failing cluster is a precise obstruction.

Neither outcome alone solves the Hubbard phase diagram or establishes a
world-level breakthrough, but either is substantially more informative than
assuming (U) is a global eliminated-sector gap.
