# Intervention and reduced-dynamics certificates: sources and derivations

This note records mathematically usable consequences of an existing positive
ground-energy certificate. It does not claim that a ground-energy certificate
alone produces a closed dynamical model.

## 1. Energy certificate to residual-observable bounds

Assume a finite-dimensional Hilbert space, a normalized state `rho`, and a
certified operator inequality

\[
 H-LI \succeq F^\dagger QF + R, \qquad Q\succeq0,\ R\succeq0,
\]

where `F` is a column of operators and `F^\dagger QF` means
`sum_ab Q_ab F_a^\dagger F_b`. If `Tr(rho H) <= L+epsilon`, positivity gives

\[
 \operatorname{Tr}(\rho F^\dagger QF)\le\epsilon.                 (1)
\]

For an operator `A=c^\dagger F` whose coefficient vector lies in the support
of `Q`, define

\[
 \kappa(c)=c^\dagger Q^+c,
\]

with `Q^+` the Moore--Penrose inverse. The block Cauchy--Schwarz inequality,
or equivalently `A^\dagger A \preceq \kappa(c)F^\dagger QF`, gives

\[
 |\operatorname{Tr}(\rho A)|\le
 \sqrt{\operatorname{Tr}(\rho A^\dagger A)}
 \le \sqrt{\kappa(c)\epsilon}.                                (2)
\]

For a matrix element, the directly safe one-sided statement is
`|<psi_1|A|psi_2>| <= ||psi_1|| sqrt(kappa(c) epsilon_2)` when the ket has
`<psi_2|F^\dagger QF|psi_2> <= epsilon_2`; a symmetric bound requires a
separate two-sided factorization or additional hypotheses. For a mixed-state
expectation, (2) is the directly safe version.

If `Q` is singular, coefficient support is an exact gate: components in
`ker(Q)` require separate certified annihilation relations or must be charged
as an unbounded residual. This is the same distinction as the project's exact
sector quotient: algebraic nullity does not imply membership in the accepting
SOS ideal.

For a residual commutator, write `[H,O]-K = c^\dagger F` and certify a bound
for `K` separately. Then (2) bounds the unresolved part of the Heisenberg
derivative. In particular, if the projected observable model uses
`dot O_red = i[H,O_red]` and the omitted term is `A`, its instantaneous
expectation error is at most `sqrt(kappa epsilon)` under the stated energy
condition. This is an observable-specific bound, not an operator-norm bound.

The derivation is elementary and does not require determinant enumeration once
the certificate and the coefficient/frame norm `kappa` have been constructed.
The expensive unresolved tasks are constructing `F,Q`, proving the operator
inequality, and controlling `epsilon(t)` after interventions.

## 2. Duhamel and residual propagation

If the retained Spectra certificate is instead in signed-residual form
`H-bI >= -eta I`, set `L=b-eta`. Then `H-LI >= 0`; the positive remainder can
be included as `R+eta I` in the factorization. The intervention bound must use
the energy excess above this shifted `L`, and must not treat a signed remainder
as positive before making this shift.

For an exact generator `G` and reduced model `G_r`, variation of constants gives

\[
 e^{tG}x(0)-x_r(t)
 = e^{tG}(x(0)-x_r(0))-
 \int_0^t e^{(t-s)G}\,r(s)\,ds,
\quad r=\dot x_r-Gx_r.
\]

Thus any verified induced norm bound on the propagator and residual gives a
finite-time certificate. For unitary Heisenberg evolution the operator norm is
norm-preserving. The adjoint of a completely positive trace-preserving map is
unital and contracts operator norm. Thus a Hermitian driven Hamiltonian does
not require a generic growth factor merely because it is time-dependent. A
verified growth factor is needed for genuinely non-contractive generators or
approximations. The project's v7
residual checker already implements this structure for piecewise polynomial
observables, including jumps and integrated residuals; this note does not
replace that checker.

For a time-dependent Hamiltonian `H(t)`, the same identity holds with the
time-ordered propagator. The energy certificate must then be uniform in time,
or the state energy must be bounded by a separate work/drive estimate. A fixed
ground-energy lower bound for `H(0)` cannot be silently reused under arbitrary
driving.

## 3. A concrete block-Krylov route

For a retained operator/state block `V`, form exact action blocks
`V, HV, H^2V, ...` and orthogonalize them in the certificate's exact metric
`G`. Retain only new quotient directions after each action and charge omitted
action as a residual. If `T_lm` is the exact Gram of action blocks and
`K_l=V^dagger H^l V`, the Schur residual Gram

\[
  R_{lm}=T_{lm}-K_l^\dagger G^+K_m
\]

is positive semidefinite when all blocks use one common exact inner product. It
is the squared norm of the action orthogonal to the retained span. This gives a
principled adaptive score and residual certificate; long-time error still
requires Duhamel propagation.

Block Krylov, Arnoldi, and Lanczos are established numerical constructions for
large eigenvalue and matrix-function problems. See Y. Saad, *Numerical Methods
for Large Eigenvalue Problems*, 2nd ed., SIAM (2011), including its chapters on
block Krylov methods and Gram matrices:
[author-hosted text](https://www-users.cse.umn.edu/~saad/eig_book_2ndEd.pdf).

The Spectra-specific opportunity is the exact rational, sector-aware metric and
Schur residual, not the Krylov idea itself. Verify `G>=0`, use `G^+` only on its
supported range, and independently replay the omitted-action norm. A small
reduced Gram or fitted factor is not a certificate unless this residual is
charged.

## 4. Audit of the exact uniform-control envelope

The implementation in `exact.py` is sound under its stated finite-dimensional
interpretation, with the following precise conditions.

Let the columns of `V` have Gram metric `G=V^dagger V>0`, let
`K_l=V^dagger H_l V`, and `A_l=G^{-1}K_l`. For Schrödinger coefficient
dynamics `dot c=-i A(u)c`, with real Hermitian `H_l` and real controls, the
diagonal entries of `A_l` contribute only phases. Therefore the componentwise
comparison matrix may set `C_ii=0` and use

`C_ij >= |A_0,ij| + sum_l a_l |A_l,ij|` for `i != j`.

The componentwise inequality `|dot c| <= C|c|` then gives
`|c(t)| <= exp(Ct)|c(0)|` for measurable controls in the declared box. This
would be false for a non-Hermitian generator or without the factor `-i`; in
those cases diagonal real parts must be included.

For each basis column `j`, the residual Gram assembled by the code is

`R_lq(j)=< (H_l V - V A_l)e_j, (H_q V - V A_q)e_j >`.

The full block matrix `[R_lq(j)]_{l,q}` is PSD. Individual cross blocks or
scalars `R_lq(j)` need not be PSD and need not be nonnegative. The control
residual square is a convex quadratic in `u`, so its maximum over a box occurs
at a vertex; enumerating all vertices is exact (though exponential in the
number of controls). Dividing by `sqrt(G_00)` correctly normalizes the initial
state `V e_0` when the reported target is the normalized state.

The Taylor integration in `integrate_envelope` is also conservative: it sums
the positive terms through the requested order and bounds the remaining
exponential tail geometrically, refusing when the ratio test is not valid.
This is an exact rational enclosure, not a numerical time integration.

The resulting `delta <= r^T integral_0^t exp(Cs)e_0 ds` bounds the state-vector
error in the `V`-induced norm only when the residual radius and propagation are
interpreted consistently. It does not certify observables outside the retained
state span. An observable conversion bound, or the observable-specific energy
bound in section 1, is required for that claim.

The previous block-residual wording should be read accordingly: the assembled
matrix is PSD; an individual nonsymmetric cross entry `R_lq` is not itself a
PSD object.

## 5. Audit of the molecular pulse trajectory checker

The trajectory theorem is sound for the declared scope: a fixed normalized
initial vector in the selected span, a specified piecewise-constant control
protocol, rational polynomial segment coefficients, and a finite-dimensional
Hermitian Hamiltonian. For each segment the checker forms the exact residual
`i V p'(theta) - duration*(H(u)-shift I)V p(theta)`, integrates its squared
Euclidean coefficient norm as a polynomial in `theta`, and takes an outward
square-root bound. Every segment interface mismatch is charged as a jump. The
phase shift is valid because subtracting a scalar multiple of the identity only
changes the global phase; it must be used consistently in both the derivative
and action terms, as the code does.

The final occupation conversion is valid for the normalized endpoint: the
orbital occupation operator has spectrum in `[0,2]`, so a normalized state-vector
error `d` changes its expectation by at most `2d`. The code separately bounds
endpoint normalization and then uses `margin=2*normalized_error`. This supports
the reported interval for that one occupation and protocol, not arbitrary
observables or controls. It also does not include physical-model error and does
not assert that the initial snapshot is the exact ground state.

An earlier inspection reported a missing `action_denominators`/
`integer_actions` integration. That reading was stale: the current
`exact_moments` return object supplies both fields, and the trajectory replay
tests and receipts have passed. The integration issue is resolved.

The construction is conventional reduced-basis/POD or Galerkin residual
certification with an unusually strict exact-rational replay and a molecular
fermionic action. Its potential contribution is the combination of that exact
certificate with Spectra's sector-aware operator representation; the pulse
result alone is not evidence of general dynamical closure, uniform-control
success, scalable discovery, or superiority to Krylov/POD methods.

## 6. What projection and memory methods contribute

Zwanzig's original projection-operator work derives an exact equation for
relevant variables with a memory convolution and an orthogonal (“random-force”)
term. It therefore supplies the correct structural template for a reduced
intervention engine: retained observables, memory kernel, and certified forcing
residual must be treated together. It does not supply a cheap kernel or a
finite error bound for a molecular SOS certificate.

Primary sources:

* R. Zwanzig, “Ensemble Method in the Theory of Irreversibility,” *J. Chem.
  Phys.* 33, 1338 (1960), DOI
  [10.1063/1.1731409](https://doi.org/10.1063/1.1731409).
* R. Zwanzig, “Memory Effects in Irreversible Thermodynamics,” *Phys. Rev.*
  124, 983 (1961), DOI
  [10.1103/PhysRev.124.983](https://doi.org/10.1103/PhysRev.124.983).
* R. Zwanzig, K. S. J. Nordholm, and W. C. Mitchell, corrected derivation,
  *Phys. Rev. A* 5, 2680 (1972), DOI
  [10.1103/PhysRevA.5.2680](https://doi.org/10.1103/PhysRevA.5.2680).

The 1972 correction matters as a warning: projection derivations can omit
fluctuation terms. Any implementation must retain the exact memory/random-force
remainder or certify its omission.

## 7. Energy-constrained norms and reduced-model scope

Energy-constrained channel norms are the appropriate language when a global
operator norm is too strong for infinite-dimensional or effectively unbounded
systems. Van Luijk proves energy-limitedness criteria for Markovian dynamics
and continuity bounds in energy-constrained norms. Shirokov develops the
energy-constrained diamond norm and its convergence/continuity properties.
These results justify an energy-budgeted validation metric, but do not turn a
finite-orbital molecular certificate into a continuum certificate.

* L. van Luijk, “Energy-limited quantum dynamics,”
  [arXiv:2405.10259](https://arxiv.org/abs/2405.10259).
* M. E. Shirokov, “Energy-constrained diamond norms and their use in quantum
  information theory,” [arXiv:1706.00361](https://arxiv.org/abs/1706.00361).

For reduced-basis methods, a residual divided by a coercivity lower bound is a
standard a posteriori pattern. A primary example develops a numerically stable
residual estimator using an orthonormal residual basis:

* Andreas Buhr, Christian Engwer, Mario Ohlberger, and Stephan Rave,
  “A numerically stable a posteriori error estimator for reduced
  basis approximations of elliptic equations,”
  [arXiv:1407.8005](https://arxiv.org/abs/1407.8005).

The analogy is useful but limited. In Spectra, `Q` and the energy certificate
play the role of a verified positive form; the exact coefficient/frame norm is
the analogue of the residual estimator constant. Unlike coercive elliptic
problems, a Hamiltonian may have degeneracies, singular forms, and time-dependent
energy injection.

## 8. Concrete implementation test

For each candidate intervention or reduced observable:

1. Export the exact rational `F,Q,R` certificate and verify
   `H-LI-F^\dagger QF-R \succeq0` independently.
2. Compute `kappa=c^\dagger Q^+c` exactly or with an outward rational enclosure.
3. Record the state energy excess `epsilon` and prove it for every time or
   intervention window used.
4. Report the observable residual bound `sqrt(kappa*epsilon)` and, for a full
   trajectory, integrate the independently bounded commutator/Duhamel residual.
5. If a drive changes `H`, regenerate or uniformly parameterize the certificate;
   do not reuse the static bound by assumption.

The first decisive experiment is therefore an observable-specific intervention
on an existing certified molecular case: compare a reduced commutator model
against exact propagation, while charging construction, `kappa`, state-energy
tracking, and the independent residual replay. Success would establish a
reusable certified observable layer; it would not yet establish closure for all
observables, arbitrary controls, or dissipative dynamics.

## 9. Audit of the robust control-neighborhood extension

The robust extension has a sound Duhamel comparison under its declared scope.
For two Hermitian Hamiltonian histories with the same initial normalized state,

\[
 \|\psi_u(t)-\psi_v(t)\|
 \le \int_0^t\|(H(u)-H(v))\psi_v(s)\|ds
 \le \sum_\ell\|H_\ell\|\int_0^t|u_\ell-v_\ell|ds .
\]

Adding this allowance to the independently replayed nominal trajectory bound
is therefore safe. The original amplitude-box check remains separate from the
integrated deviation budget.

For an even, number-conserving CAR polynomial whose support has at most eight
modes, permuting the supported modes gives a local Fock matrix tensored with
identity on the complement. Its maximum absolute row sum bounds the operator
2-norm for the Hermitian local matrix. For larger support, the coefficient
one-norm is a valid triangle bound because each normalized CAR monomial has
operator norm at most one. These are conservative support-based bounds; they do
not establish cheap construction for large support.

The all-measurable-perturbation quantifier is valid only for the stated
integrated absolute budgets and common initial state. It does not cover a
changed initial state, perturbations outside the original amplitude box, or
perturbations to terms not included in the Hamiltonian decomposition.

The current trajectory checker strengthens this to a uniform-in-time bound
without needing explicit prefix records. Every prefix residual integral is at
most the full-segment `L1` residual bound `sqrt(I_k)` by Cauchy--Schwarz on a
subinterval of `[0,1]`, and every prefix contains only a subset of the charged
nonnegative jumps and segment contributions. Therefore the single total
`sum_k (jump_k+sqrt(I_k))` dominates every prefix. The implementation now also
computes an exact Bernstein enclosure of
`||Vp(theta)||^2/G00` on every segment, including both sides of jumps. Taking
the worst deviation of the enclosed norm from one and adding it to the total
phase-aligned bound gives a valid uniform normalized state-error bound, provided
the enclosure's Bernstein degree and outward rounding are accepted. Prefix
records and partial residual integrals would only tighten the result; they are
not required for soundness of this conservative uniform certificate.
## 10. Scope of the fixed-linear-subspace obstruction

The new `MATHEMATICS.md` obstruction is correctly limited to invariance under
*all* number-conserving one-body spin-orbital controls in a fixed `N` sector.
It does not apply to the two tested molecular controls, to spin-invariant-only
control families, or to adaptive/nonlinear representations. The occupation
projector plus hopping argument establishes irreducibility of that larger
control algebra, but cannot be used to infer that a particular two-control
family has no compact representation. Nor does it rule out query-specific,
trajectory-specific, nonlinear, tensor-network, or adaptively enriched models.

## 11. Audit of the compiled joint-action kernel

`kernel.py` constructs the joint integer action matrix from the original
Hamiltonian, proposal columns, and all reached labels before discarding the
configuration data. Its Gram is exactly

\[
 J=(V,H_0V,DV,WV)^\dagger(V,H_0V,DV,WV),
\]

so all mixed action terms are retained. Separate action denominators are folded
into a common query denominator with `lcm`; the duration and phase-shift factors
in `residual_integral` then cancel against the final denominator as required.
The compiled path uses the same metric for jumps and endpoint norms, and the
observable Gram for the occupation quotient. Matching direct and compiled H6
receipts, plus the nonorthogonal/control/shift tests, supports equivalence of
the replay paths.

This is a reusable compiled verifier: after the charged exact construction,
later queries use only small exact matrices and retain zero configuration
records. It is not enumeration-free construction, since kernel creation visits
all supplied/reached configurations and stores a dense joint Gram. The small
Gram is not self-authenticating; trust still depends on fresh original-input
construction and fixture/proposal/orbital bindings. Any amortized speedup
therefore applies to many queries sharing one compiled problem and proposal,
not to a cheap one-shot many-body solver or a universally small proof.
