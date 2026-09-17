# Certified molecular reduction under declared interventions

This is a concrete finite-model investigation of intervention-valid reduction.
The mathematics uses established Galerkin, residual, comparison-system, and
snapshot-reduction ideas. Neither the identities nor the numerical results are
claimed to solve general many-body computation or to be a new field-wide method.

## Domain and initial state

Work in a fixed-particle sector of a finite fermion Fock space. The Hamiltonian
and controls are real Hermitian CAR polynomials supplied as exact rational
coefficients:

\[
 H(u,t)=H_0+\sum_{\ell=1}^q u_\ell(t)H_\ell,
 \qquad |u_\ell(t)|\le a_\ell.
\]

Atomic units set hbar=1. For the current molecular tests, q=2. On adjacent
spatial orbitals p,q, the controls are

\[
 D=n_{p\alpha}+n_{p\beta}-n_{q\alpha}-n_{q\beta},\qquad
 W=\sum_{\sigma}(a^\dagger_{p\sigma}a_{q\sigma}
                +a^\dagger_{q\sigma}a_{p\sigma}).
\]

These are noncommuting one-body perturbations to the existing interacting
molecular Hamiltonian. They are diagnostic Hamiltonian controls; no claim is
made that an experimental electric field implements these precise operators.

Let V contain real rational state columns supported on supplied valid
configurations. The declared normalized initial state is

\[
 \psi_0=Ve_0/\sqrt{G_{00}},\qquad G=V^\dagger V>0.
\]

It is not silently identified with the unknown exact molecular ground state.
An energy upper bound does not supply a state-distance bound without further
information. No laboratory preparation fidelity or physical-model error is
included. Number conservation is verified. Total-spin purity is not presumed.

The current implementation stores these columns on configurations. It may
visit or retain the entire magnetic sector. A low reduced dimension is therefore
not evidence of enumeration-free construction. Counts of both retained and
reached labels are part of the acceptance receipt.

## Exact projected model and omitted action

Define

\[
 K_\ell=V^\dagger H_\ell V,\quad A_\ell=G^{-1}K_\ell,
 \quad B_\ell=H_\ell V-VA_\ell.
\]

The exact reduced equation is

\[
 i\dot c=A(u)c,\quad c(0)=e_0,
 \qquad \widetilde\psi=Vc/\sqrt{G_{00}}.
\]

Although A need not be symmetric in nonorthogonal coordinates,
A^dagger G=GA. Thus c^dagger Gc is preserved and the lifted reduced state
is normalized. Every scalar diagonal entry of A is real in this implementation.

The mixed residual matrices obey

\[
 B_\ell^\dagger B_m
 =V^\dagger H_\ell H_mV-K_\ell G^{-1}K_m.
\]

The assembled block matrix is positive semidefinite. An individual cross block
need not be Hermitian or positive. Mixed control terms are retained exactly;
replacing the residual square by a sum of unmixed squares can underestimate it.

For any control history, unitary variation of constants gives

\[
 \|\psi(t)-\widetilde\psi(t)\|
 \le\frac1{\sqrt{G_{00}}}\int_0^t\|B(u(s))c(s)\|\,ds.
\]

No excitation gap, adiabatic assumption, or ground-state uniqueness is needed.
The omitted action and its construction are still essential costs.

## A uniform all-waveform bound

Set C_ii=0 and, for i != j,

\[
 C_{ij}=|(A_0)_{ij}|+\sum_\ell a_\ell|(A_\ell)_{ij}|.
\]

The real diagonal part generates phases and contributes no growth to |c_i|.
The upper Dini derivative satisfies D^+|c| <= C|c|, including at zeros.
The positive comparison principle then implies

\[
 |c(t)|\le e^{Ct}e_0.
\]

Each column's residual norm is convex in u. Therefore its maximum on the
control box is achieved at a vertex. Exact arithmetic at all vertices gives
r_j >= max_u ||B(u)e_j||/sqrt(G_00). Consequently

\[
 \boxed{\|\psi(t)-\widetilde\psi(t)\|
 \le r^T\int_0^t e^{Cs}e_0\,ds.}
\]

This covers arbitrary bounded measurable controls in the box, including rapid
switching. No finite sampling of control histories is used to justify that
quantifier. Coefficients of C and r are rounded outward on a rational grid.

For x=T||C||_1 and truncation order p with x/(p+3)<1, the omitted positive
Taylor tail of the integral is at most

\[
 r_{\max}T\,
 \frac{x^{p+1}}{(p+2)!}\frac1{1-x/(p+3)}.
\]

The exact retained sum and this tail prove the printed bound. Initial tests
showed this uniform bound was too weak at the frozen molecular target. That is
a failure of those certificates, not a proof that no small valid model exists.

## A joint trajectory residual that preserves interference

For a specified piecewise-constant protocol, construct a rational complex
polynomial p(theta) on each segment, theta in [0,1], of duration h. A real
constant energy shift s removes a global phase. Let

\[
 v(t)=Vp(\theta)/\sqrt{G_{00}},\quad
 r(\theta)=iVp'(\theta)-h(H(u)-sI)Vp(\theta).
\]

The accepting checker reconstructs every coefficient of this vector from
the original CAR Hamiltonian, controls, and integer state columns. It does not
trust the projected numerical generator or its eigensolver. In particular,
derivative terms, control terms, and basis columns are combined before taking
the squared norm.

Writing r(theta)=sum_k r_k theta^k, the integrated squared residual is the
exact nonnegative rational

\[
 I=\frac1{G_{00}}\sum_{k,j}
   \frac{\operatorname{Re}(r_k^\dagger r_j)}{k+j+1}.
\]

By Cauchy--Schwarz on an interval of length one,

\[
 \int_0^1 \|r(\theta)\|/\sqrt{G_{00}}\,d\theta\le\sqrt I.
\]

The segment duration is already inside r. Multiplying by h again would be
incorrect. The implementation uses a common integer denominator for the
operator, duration, polynomial, and input-state normalization, then an outward
rational square root.

All initial and inter-segment jumps are reconstructed and charged. If J_k is
the norm of jump k, the endpoint phase-aligned state error is bounded by

\[
 d=\sum_k J_k+\sum_k\sqrt{I_k}.
\]

This includes reduced-model error, numerical time propagation, rational
rounding, and omitted configurations. It needs no numerical integrator to be
trusted. A successful replay means the inequality was proved; meeting the
requested numerical tolerance is a separate field in the receipt.

The polynomial endpoint need not have norm one. Its norm is evaluated exactly
and enclosed outward. Normalizing it adds at most |1-||v(T)||| to d. This yields
the printed normalized state-distance bound d_norm.

For any Hermitian observable whose spectrum lies in [a,b], the resulting
expectation error is at most

\[
 (b-a)\min(1,d_{\rm norm}).
\]

This follows from the pure-state trace-distance bound. The implemented example
is spatial occupation, whose spectrum lies in [0,2]. Its reduced expectation is
an exact rational quotient. The interval is widened by 2*d_norm and intersected
with [0,2]. Subtracting the exactly known initial expectation can establish a
signed change. The state certificate can support other bounded observables,
but each new observable's matrix elements and spectral range must be checked.

## Uniform time and control neighborhoods

The total sum of full-segment residual bounds and all jump norms also bounds
every prefix of a trajectory. For a partial segment ending at theta <= 1,
the integral of the residual norm is at most sqrt(theta) times the square root
of its partial squared integral, hence at most sqrt(I_k) for the whole segment.
All terms in the accumulated bound are nonnegative. Detailed prefix records
could tighten this estimate but are not necessary for its validity.

On each segment compute the exact polynomial

\[
 f(\theta)=\frac{p(\theta)^\dagger V^\dagger Vp(\theta)}{G_{00}}
 =\sum_{j=0}^{d}a_j\theta^j.
\]

Its Bernstein coefficients on [0,1] are

\[
 b_k=\sum_{j=0}^{k}a_j\frac{\binom{k}{j}}{\binom{d}{j}}.
\]

The Bernstein basis is a nonnegative partition of unity, so
min(b_k) <= f(theta) <= max(b_k). A strictly positive lower bound establishes
that the polynomial can be normalized throughout the segment. Outward square
roots and the maximum normalization correction across every segment yield a
uniform-in-time normalized-state error. If the sufficient positivity check
fails, the implementation withholds this uniform normalized claim even when
an endpoint certificate exists.

Let u be a checked nominal protocol and v any measurable perturbed protocol.
Both must satisfy the originally declared amplitude box. Unitarity gives

\[
 \|\psi_v(t)-\psi_u(t)\|
 \le\sum_\ell\|H_\ell\|\int_0^t|v_\ell(s)-u_\ell(s)|ds.
\]

Thus integrated absolute deviation budgets rho_l enlarge either the endpoint
or the uniform nominal state bound by sum_l ||H_l|| rho_l. They cover a
continuous family of waveforms, not just sampled perturbations. For a uniform
amplitude deviation delta_l over horizon T, rho_l=T delta_l. For the two tested
controls the exact full-Fock norm bounds are both 2.

The implementation proves these norm bounds on the local CAR support. An even
operator is unitarily equivalent after a mode permutation to its local Fock
matrix tensor identity. Its maximum absolute row sum bounds the Hermitian
spectral norm. Only 16 local labels are used for each four-mode control. A CAR
coefficient one-norm fallback avoids local enumeration for larger supports.
This small local calculation does not remove the separate many-body dependency
of the dynamical reduction.

## Reusable exact action kernel

Define the joint column operator and its Gram matrix

\[
 Z=(V,H_0V,H_1V,\ldots,H_qV),\qquad J=Z^\dagger Z.
\]

Every coefficient of a polynomial trajectory residual is Z times a small
coefficient vector. Therefore all mixed residual inner products are exact
quadratic or bilinear forms in J. The ordinary metric and the selected
observable matrix provide the jump, normalization, and expectation terms.
For q=2 the joint dimension is 4r. No commutation assumption is used; all
interference terms are retained.

The implementation stores integer column numerators with separate exact
denominators. Its stored integer Gram is converted to the displayed rational J
by the corresponding diagonal denominator factors. Those factors enter every
query; raw integer action columns are not mistaken for normalized physical ones.

`CheckedKernel` constructs J from original CAR actions with integers. It
retains no configuration records for later queries. The query API uses only
the verified small matrices, bindings, and local controls. It does not load an
unverified matrix file as an accepting shortcut. The bundle path discards the
original configuration arrays before querying.

Construction still visits every supplied and reached configuration. The
kernel is compiled verifier state, not a standalone small proof that J equals
its many-body definition. Its initial verification and the proposal's discovery
must be charged. The value demonstrated here is reuse of already checked
quadratic data, not enumeration-free discovery.

## Numerical proposals and finite inverse control

Neither Taylor propagation nor matrix-exponential interpolation is trusted by
acceptance. Stable Taylor adapts the step to the numerical projected generator.
Chebyshev interpolation samples the reduced matrix exponential, rounds the
Chebyshev coefficients to integers, then changes to power coefficients using
the exact recurrence for T_k(2 theta - 1). The latter avoids unstable floating
conversion followed by pretending that cancellation is exact. Every resulting
polynomial is subjected to the same independent rational residual test.

The inverse experiment evaluates a finite menu of 729 three-stage protocols in
the reduced model, chooses the largest predicted population increase, and
checks that candidate in the original Hamiltonian. A target is accepted only
when the lower endpoint of the rigorous change interval exceeds it. Numerical
ranking is not a proof of a global optimum. The delivered policy is a finite-
Hamiltonian control sequence, not a molecular geometry or laboratory synthesis
policy. The initial state remains the specified mathematical vector.

## Discovery and adaptive enrichment

The first selectors kept low-energy eigenstates, then states strongly coupled
to the initial state by the controls. These did not meet the uniform target.
The next selector uses standard real snapshot/POD directions from four constant
corner controls. It inherits MPS-guided configuration selection and uses dense
selected-configuration eigensystems. These are charged discovery dependencies.

Residual enrichment constructs the actual vectors

\[
 [H(u)V-VA(u)]\,e^{-itA(u)}e_0
\]

at the frozen training corners and times. Leading independent real directions
are added to V, followed by a new exact trajectory check. This is ordinary
residual-driven reduced-basis enrichment; its value here is measured accepted
error improvement. The enlarged constructor applies sparse CAR actions and
does not diagonalize an enlarged configuration matrix. It does materialize all
reached vector labels. The initial column is preserved exactly.

The protocol-specific certificate has a narrower quantifier than the uniform
control-box theorem. A passing pulse must never be presented as a proof for all
waveforms. Training trajectories, evaluation pulses, and adaptive retests are
identified in the protocol and receipt files.

## A precise obstruction to one overly strong interpretation

No nonzero proper fixed linear state subspace is invariant under *all*
number-conserving one-body controls on m spin orbitals in an N-particle sector.

An elementary proof is useful. Invariance under n_i implies invariance under
products of n_i and (I-n_i), including each occupation projector. Applying a
projector to a nonzero coefficient of any vector in the invariant subspace
puts at least one occupation basis state in that subspace. Hermitian hopping
controls a_p^dagger a_q+a_q^dagger a_p move an occupied q to an empty p.
Repeating such moves reaches every N-particle configuration. The subspace must
therefore have dimension binomial(m,N).

This is a standard irreducibility fact with the proof supplied here, not a new
hardness theorem. Its assumptions are stronger than the two controls tested
above; it is not an obstruction for those particular molecular cases. It also
does not exclude nonlinear or adaptive representations: Slater determinants,
for example, remain compact under purely quadratic dynamics even though their
linear span fills the particle sector. Thus a universal *fixed linear closure*
would be the wrong success requirement. Approximate query-specific and adaptive
representations remain open possibilities.

## Remaining constructive gap

These identities provide an accepting boundary for intervention models. They
do not discover a uniformly cheap many-body representation, remove the global
configuration dependency of this prototype, control continuum/model error,
establish finite-temperature behavior, or construct a physical synthesis route.
The residual-driven directions are evidence only where exact acceptance and
complete accounting demonstrate a gain. See SOURCES.md for the primary-source
context and mathematical review.
