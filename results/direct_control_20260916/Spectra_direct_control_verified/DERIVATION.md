# Direct collective-control certificate and its limits

The accepted result is a strong-control component. The original low-amplitude
H8 problem and the general many-body problem remain unresolved. Its accepting
calculation uses integer/rational arithmetic. The separate tensor diagnostics
have additional floating-arithmetic assumptions.

## Model and control algebra

For spatial orbitals with even spin indices p and q, define

\[
D=n_p+n_{p+1}-n_q-n_{q+1},\qquad
W=\sum_{\sigma=0}^1(a_{p+\sigma}^\dagger a_{q+\sigma}
+a_{q+\sigma}^\dagger a_{p+\sigma}),
\]

\[
A=\sum_{\sigma=0}^1(a_{q+\sigma}^\dagger a_{p+\sigma}
-a_{p+\sigma}^\dagger a_{q+\sigma}),\qquad Y=iA.
\]

Exact CAR identities give

\[
[W,D]=2A,\qquad [W,A]=2D,\qquad A^\dagger=-A.
\]

Both D and W have norm at most two and preserve each spin population. The
checker verifies the polynomial identities. Real initial tensors imply
\(\langle Y\rangle_0=0\). On H8, `(p,q)=(6,10)` and
\(d_0=\langle D\rangle_0=1.390133812147592\ldots\).

For nominal \(H=H_0+vW\), `u=0`, `v>0`, the control-only observable is

\[
D_c(t)=e^{ivWt}De^{-ivWt}=D\cos(2vt)+Y\sin(2vt).
\]

Thus its initial-state expectation is \(d_0\cos(2vt)\).

## Collective treatment of the interacting drift

Let \(\mathcal L(X)=i[H,X]\) and \(\mathcal L_c(X)=i[vW,X]\). Duhamel gives

\[
e^{T\mathcal L}D-e^{T\mathcal L_c}D
=\int_0^T e^{(T-s)\mathcal L}i[H_0,e^{s\mathcal L_c}D]ds.
\]

Unitary conjugation preserves norm. If
\(C_D\ge\|[H_0,D]\|\), \(C_Y\ge\|[H_0,A]\|=\|[H_0,Y]\|\), then

\[
|\langle D(T)\rangle-d_0\cos\theta|
\le C_D\int_0^T|\cos(2vt)|dt+C_Y\int_0^T|\sin(2vt)|dt,
\quad \theta=2vT.
\]

For \(\pi/2<\theta<\pi\), this is

\[
\boxed{\epsilon_{\rm drift}=
\frac{C_D(2-\sin\theta)+C_Y(1-\cos\theta)}{2v}.}
\]

All original interactions remain present through the two commutator bounds.
No trajectory, embedding, or full sector matrix is needed. This succeeds when
the control acts before drift accumulates appreciably; it does not resolve
arbitrary strongly correlated dynamics.

## Exact computational inputs

Convert the CAR polynomials exactly to products of `I, X, J=-iY, Z`. Every
product has norm one. Combine identical products before summing absolute
rational coefficients. A quartic CAR word generates at most 16 choices before
combination. This gives norm upper bounds without configuration enumeration.

For H8:

\[
C_D=19.128975101863\;\mathrm{Ha},\qquad
C_Y=18.676738267037\;\mathrm{Ha},
\]

with 2,236 and 1,992 nonzero products respectively. These are conservative
bounds, not measured spectral norms.

The initial expectation is an integer MPS contraction divided by the exact
integer norm. The common tensor denominator cancels. Charge flow, positive
norm, and fixture/state hashes are checked. H8 needs 962,434 integer
multiplications and at most 2,098 simultaneous transfer entries. Its supplied
MPS has 10,093 nonzero tensor entries and maximum bond 144. Cost still depends
on those bonds and integer sizes; their initial discovery is inherited.

For rational \(\theta=25/8\), alternating Taylor sums enclose sine and cosine.
An exact Machin-formula enclosure of pi verifies the required angle range.
Floating endpoint displays are not used in acceptance.

## Simultaneous uncertainty budgets

Pointwise errors at most a in both control coefficients give
\(\|\delta uD+\delta vW\|\le4a\). Trace-norm stability therefore contributes
at most \(16Ta\) to the observable. Initial trace distance r0, using half the
trace norm, contributes \(4r_0\).

The declared noise model is

\[
\dot\rho=-i[H(t),\rho]+\sum_j\gamma_j(t)(Z_j\rho Z_j-\rho),
\quad Z_j=I-2n_j,\quad\gamma_j\ge0,
\quad\int_0^T\sum_j\gamma_jdt\le\Gamma.
\]

Since \(\|Z_j\rho Z_j-\rho\|_1\le2\), contractivity gives allowance
\(4\Gamma\). Altogether,

\[
\boxed{\epsilon_{\rm robust}=16Ta+4r_0+4\Gamma.}
\]

This covers every continuous-time profile satisfying these assumptions, not
only sampled profiles. Intersect the rational enclosure of
\(d_0\cos\theta\pm(\epsilon_{\rm drift}+\epsilon_{\rm robust})\)
with `[-2,2]`.

Using the source short-case budgets `a=1/1000 Ha`, `r0=1/2000`, and
`Gamma=1/1000`, doubling v from 1/2 selects

\[
v=64\;\mathrm{Ha},\quad T=25/1024,
\qquad \langle D(T)\rangle\le-0.7953372198008\ldots<-0.6.
\]

The other endpoint is approximately -1.9845476869817. The receipt contains
exact rational endpoints. The amplitude is 128 times the original cap. This
changes a material task constraint and is explicitly a separate regime.

## Exact limitation of the weak-control envelope

This section bounds the **upper endpoint of the proof rule**, not the actual
observable. Apply `[H0,D]` and `[H0,A]` to one balanced-spin determinant,
integer label 255. Their direct CAR images have support 123 and 134. The exact
squared image norms give rational lower bounds

\[
L_D>0.6188411437470,\qquad L_Y>0.7687226806648.
\]

Any valid commutator norm bound must exceed these witnesses. With
\(v_{\max}=1/2\), set \(a=L_D/(2v_{\max})\), \(b=L_Y/(2v_{\max})\).
For every \(0<|v|\le v_{\max}\), use \(\theta=2|v|T\). Dropping nonnegative
robustness allowances, the envelope's upper endpoint obeys:

* If \(0\le\theta\le\pi/2\), it is at least zero, since the ideal expectation
  and drift allowance are nonnegative.
* If \(\pi/2\le\theta\le\pi\), it is at least
  \(2a+b-\sqrt{(d_0-b)^2+a^2}>1.1294\). This minimizes a linear form in sine
  and cosine over the entire unit circle, a relaxation of the actual quadrant.
* If \(\theta\ge\pi\), accumulated absolute sine and cosine integrals are
  each at least two. The endpoint is at least \(2a+2b-d_0>1.3849\).

Zero amplitude also cannot make this rule give a negative upper bound. Capping
the endpoint at two leaves it nonnegative. Thus **no duration** at the original
amplitude cap can certify the negative target through this rule, even with
exact commutator norms. At the fixed angle 25/8, including the uncertainty
budget gives the stronger endpoint floor 1.4308116314.

The obstruction recomputes the bound MPS's initial moment and checks the
real-tensor/zero-Y premise. It concerns constant W-only control with `u=0` and
this norm-sum envelope. It does not constrain other proof mechanisms, physical
reachability, or the source four-phase control family.

## State-compression diagnostics

For Schmidt coefficients \(\sigma_j\) at a cut of a normalized state, every
rank-r approximation has vector-norm error at least

\[
\left(\sum_{j>r}\sigma_j^2\right)^{1/2}.
\]

Tensor canonicalization gives a small center matrix. The diagnostics bound
actual reconstruction residuals and isometry defects. Subtract their total
allowance from a checked singular-tail lower bound, and divide by an upper
bound on the source norm. Final square-root rounding is rational, but the
BLAS/IEEE floating error model remains an explicit assumption.

Canonical-order H8 has state-error lower bounds 0.0566058 at bond 64 and
0.0111295 at bond 128. In the declared local basis, the bond-64 lower bound is
0.0052285. These results concern state norm in those bases. They do not prove
that a compact representation of one observable is impossible, or that a
larger bond suffices for accurate time evolution.

## A correct possible adjoint refinement

This is recorded to prevent an invalid shortcut; it produced no accepted
weak-control result here. For a continuous candidate q, let
\(r=i\dot q-Hq\), \(e_0=\psi(0)-q(0)\), and
\(E=\|e_0\|+\int_0^T\|r\|dt\). Then

\[
e(T)=U(T,0)e_0+i\int_0^TU(T,s)r(s)ds.
\]

For Hermitian O set \(z(s)=U(s,T)Oq(T)\). The exact identity is

\[
\langle\psi(T),O\psi(T)\rangle-\langle q(T),Oq(T)\rangle
=2\operatorname{Re}\left(\langle z(0),e_0\rangle
+i\int_0^T\langle z(s),r(s)\rangle ds\right)
+\langle e(T),Oe(T)\rangle.
\]

For an approximate backward vector p, define \(s_p=i\dot p-Hp\),
\(\delta_T=\|Oq(T)-p(T)\|\), and
\(\eta(t)=\delta_T+\int_t^T\|s_p(s)\|ds\). Replacing z by p requires
remainder allowance

\[
2\eta(0)\|e_0\|+2\int_0^T\eta(t)\|r(t)\|dt+\|O\|E^2.
\]

Trajectory jumps need additional terms. Omitting this remainder is invalid.
For a forward density matrix obeying `rho_dot=-i[H,rho]`, the backward
observable's defect is `O_dot+i[H,O]`. Neither notation eliminates the need for
tight bounds on omitted errors.

## Scope and attribution

CAR/Jordan–Wigner algebra, MPS contractions, Duhamel stability and Schmidt-tail
bounds are established ingredients. This work supplies checked instances and
reusable code, not a priority claim for those principles. Initial Hamiltonian
and MPS discovery, and the local-orbital proposal, are inherited dependencies.
No result here establishes economical arbitrary many-body evolution, a
laboratory actuator, or a synthesis policy.
