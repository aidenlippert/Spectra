# Local coercive response theorem

This is a concrete route for an interacting Hubbard-type Hamiltonian. It
replaces a global scalar denominator gap by local positive pieces plus a
certified perturbation. It is a sufficient construction, not a claim that
every many-body Hamiltonian has the required decomposition.

Let (Q) be the eliminated sector and write

\[
 A=Q(H-b)Q=D+V,
\qquad D=\sum_{C\in\mathcal C}D_C.
\]

Assume each local (D_C\succeq0), with a local excitation coercivity

\[
 D\succeq d I_Q,
\qquad d>0,
\]

and certify the coupling remainder by a local overlap coloring:

\[
 \|V\|\le\kappa<d,
 \qquad\kappa\le\sum_{r=1}^{\chi}
 \sum_{C\in\mathcal C_r}\|V_C\|.
\]

The sum is the safe triangle bound. Disjoint supports do **not** make the
operator norm a maximum: for example, \\(\|\sum_i Z_i\|=N\\) although each
\\(\|Z_i\|=1\\). Any sharper extensive or boundary estimate requires an
additional commuting-spectrum or relative-form argument and must be proved
separately. For a Hubbard ladder, (D_C) may
contain (U n_{i\uparrow}n_{i\downarrow}), local chemical-potential terms,
and positive boundary penalties; (V_C) contains the hopping and residual
terms assigned to that cluster. The displayed bound is checkable from local
CAR word norms and a finite coloring, without listing determinants.

## Certified denominator bounds

The local decomposition implies

\[
 (d-\kappa)I\preceq A\preceq(d_D+\kappa)I,
 \qquad d_D=\|D\|.
\]

More usefully, if (D^{-1}) is available as a local response action and

\[
 \|D^{-1/2}VD^{-1/2}\|\le\rho<1,
\]

then the relative sandwich

\[
 (1-\rho)D\preceq A\preceq(1+\rho)D
\]

is certified, and therefore

\[
 \frac1{1+\rho}D^{-1}\preceq A^{-1}preceq
 \frac1{1-ho}D^{-1}.
\]

No square root of (A) is needed by the checker. A conservative sufficient
condition is (|V\|/d\leho). A sharper local SOS proof may certify the
relative bound directly on each overlap color.

## Local response and residual certificate

Choose a finite-depth Neumann response, **provided an independently
constructible action of (D^{-1}) exists**:

\[
 X_k=\sum_{j=0}^{k-1}(-D^{-1}V)^jD^{-1}B,
 \qquad B=QHP.
\]

Its residual is exactly

\[
 R_k=B-AX_k=(-1)^k(VD^{-1})^kB.
\]

Consequently, with (G=D^{-1/2}B),

\[
 R_k^*D^{-1}R_k
 \preceq \rho^{2k}G^*G.
\]

Using the denominator sandwich gives the fully checked response enclosure

\[
 2\operatorname{Re}(X_k^*B)-X_k^*AX_k
 +\frac1{1+\rho}R_k^*D^{-1}R_k
 \preceq B^*A^{-1}B
 \preceq
 2\operatorname{Re}(X_k^*B)-X_k^*AX_k
 +\frac1{1-\rho}R_k^*D^{-1}R_k.
\]

The finite Neumann word count alone does not make (D^{-1}) local: even a sum
of local terms can have a globally supported inverse. The acceptance checker
needs an explicit representation or certificate for every (D^{-1}) action,
in addition to the exact residual identity,
the local PSD/coercivity certificates, and the final retained inequality.
It does not need a determinant basis.

## The real remaining obligation

The theorem controls the eliminated response. It does **not** certify that the
retained Schur operator is positive. That still requires a global retained
certificate, or another local coercive decomposition applied recursively. A
claimed breakthrough must demonstrate that the retained remainder admits the
same colored local structure and that (ho<1) does not deteriorate with
system length or elimination depth.

## Adversarial limits

Local coefficient size alone is insufficient. A rank-one long-range coupling
can have small coefficients but operator norm proportional to system size;
the coloring bound must therefore count overlap and support, not coefficient
size alone. Also, (U>0) does not imply a gap on a sector containing states
with no doublons. The eliminated sector must be selected so that the local
coercivity (d>0) is actually true. Finally, if (ho\ge1), the Neumann
construction is invalid even though (A) may remain positive; this is a
diagnostic failure of the chosen local preconditioner, not a theorem that the
physics is hard.

The new mathematical promise is conditional and measurable: if local
coercivity, finite overlap coloring, and a uniform relative residual hold,
response width decays geometrically as (ho^{2k}) while construction cost
scales with local supports and depth. Establishing those conditions for a
nontrivial family of strongly correlated ladders, and closing the retained
positivity proof, remains open.
