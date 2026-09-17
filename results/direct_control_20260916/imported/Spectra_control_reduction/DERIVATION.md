# Certified reduced dynamics and robust molecular steering

## Status

This package adds a standard-library, rational acceptance path for driven quantum dynamics on the specified H8 Hamiltonian. It proves finite-time orbital-population reversal starting from the supplied rational MPS, including explicitly bounded waveform deviations, initial-state uncertainty, and phase-flip dephasing.

It is **not** a universal many-body solution or a demonstrated economical discovery method. Discovery and acceptance enumerate the 4,900-dimensional balanced-spin sector. Numerical full-model snapshots discover the reduced bases. There is no laboratory control mapping, state-preparation demonstration, exact-ground initialization claim, physical-error certificate, H12 result, or claimed mathematical priority.

A-posteriori variational dynamics bounds and offline/online reduced-basis methods have established predecessors; see SOURCES.md. The new work is the actual exact implementation, molecular certificates, robust targets, and documented failures.

## 1. Actual dynamical problem

The input is the unchanged rational electronic H8 Hamiltonian, eight electrons in sixteen spin orbitals, in canonical RHF orbitals for a 1.4 Angstrom STO-3G chain. Alpha and beta occupations are both four. The witness is the supplied charge-flow MPS with its actual integer tensors. It is not asserted to be the exact ground state.

Define

\[
D=n_6+n_7-n_{10}-n_{11},\qquad
W=a_6^\dagger a_{10}+a_{10}^\dagger a_6+
 a_7^\dagger a_{11}+a_{11}^\dagger a_7.
\]

The full model obeys

\[
i\partial_t\Phi(t)=K(t)\Phi(t),\quad
K(t)=H_0+(1851/200)I+u(t)D+v(t)W.
\]

The scalar shift only changes a global phase. Observable predictions therefore apply to H0+uD+vW. Time is in atomic units (hbar/Hartree); controls have Hartree units. The defined controls are Hamiltonian coefficients, not demonstrated electric-field, laser, geometry, or synthesis settings.

D and W are Hermitian, preserve both spin populations, do not commute, and have operator norms at most two. D's bound is immediate from its four local occupation bits. W is two commuting spin-resolved hopping operators, each of norm one; the independent local occupation oracle also verifies its Hermiticity and absolute row-sum bound.

All original H0 terms remain present. The new calculation uses no spin averaging or additional spin approximation: literal action on the invariant balanced-spin sector is checked.

## 2. Defect propagation without an exponential stability factor

Let q(t) be any piecewise continuously differentiable proposed vector path. Define

\[
r(t)=i\dot q(t)-K(t)q(t).
\]

For the exact solution with initial vector Phi0, unitarity and variation of constants give

\[
\|\Phi(t)-q(t)\|\le
\|\Phi_0-q(0)\|+
\int_0^t\|r(s)\|ds+
\sum_{\text{jumps before }t}\|q(s^+)-q(s^-)\|.
\]

This does not require q to solve a variational equation, preserve norm, be differentiable at control switches, or have a small pointwise energy error. It uses Hermiticity of the actual generator. There is no exp(||H||t) Gronwall factor. That is established unitary-stability mathematics, not a claim of a new general stability theorem.

This bound also explains the exact limitation of ground-state-only reduction. In a three-level example H0=diag(0,1,3), the subspace span{|0>,|1>} preserves the ground energy exactly. A resonant control g(exp(-3it)|2><0|+exp(3it)|0><2|) transfers |0> completely to |2> at gt=pi/2, whereas its projected control is zero. Thus exact equilibrium reduction alone does not certify interventions. The omitted control action must be bounded too.

## 3. Reduced representation and exact small defect matrices

Write q(t)=J z(t), where J is a rectangular rational matrix. It need not be exactly isometric. For each nominal control center, a real symmetric rational model h is proposed. Define

\[
E=KJ-Jh,\quad G=J^\dagger J,\quad F=E^\dagger E.
\]

Then

\[
r=J(i\dot z-hz)-Ez,
\]

so

\[
\|r\|\le\|J\|\,\|i\dot z-hz\|+\sqrt{z^\dagger Fz}.
\]

The checker constructs J's actual Gram matrix and E's actual Gram matrix, using all molecular terms. It does not accept floating projected residuals, POD singular-value tails, or a claimed Galerkin identity. The rounded h can be imperfect: its mismatch is explicitly in E.

A Gershgorin/absolute-row-sum bound on G supplies a rational enclosure of ||J||. The matrices controlling trajectory error are only r-by-r after this expensive preparation.

**The cost boundary is explicit:** this implementation computes these matrices from dense basis columns on the 4,900-state sector, not from cheap tensor contractions. For long64 the embedding alone contains 313,600 integer coefficients. The abstract formulas permit other backends, but no enumeration-free backend was obtained here.

## 4. Continuous time, not quadrature acceptance

Each time segment of length dt uses the literal polynomial

\[
z(s)=\sum_{k=0}^p c_k s^k,\quad s\in[0,1],\qquad t=t_0+dt\,s.
\]

Here p=18 in the saved examples; all coefficients have rational real and imaginary parts. They were proposed by numerical reduced dynamics, but every mismatch is recomputed.

For the representation part of the defect,

\[
I=\int_0^1 z(s)^\dagger Fz(s)\,ds
=\sum_{k,l}\frac{\operatorname{Re}(c_k^\dagger F c_l)}{k+l+1}.
\]

Since this is the integral of a squared norm,

\[
\int_{t_0}^{t_0+dt}\|Ez\|dt\le dt\sqrt{I}.
\]

The exact Gram entries are enclosed on a dyadic grid. Let Flo be the entrywise floor with grid q=2^90, so |Fij-Floij/q|<1/q. For integer polynomial coefficients with common denominator zden, the checker evaluates the integral for Flo exactly and adds

\[
\frac{\left(\sum_{k,i}(|\Re c^{\rm int}_{k,i}|+|\Im c^{\rm int}_{k,i}|)\right)^2}
{q\,zden^2}.
\]

This is conservative because it bounds the integral of ||z(s)||_1^2 on 0<=s<=1. Imaginary mixed terms cancel for the real symmetric F. An exact lcm of 1,...,2p+1 evaluates the remaining power integrals. No sample grid is the acceptance argument.

For the reduced integration defect, write

\[
i\dot z-hz=\sum_k d_k s^k.
\]

Literal rational coefficient reconstruction and the triangle inequality give

\[
\int\|J(i\dot z-hz)\|dt
\le dt\,\|J\|\sum_k\frac{\|d_k\|}{k+1}.
\]

The final coefficient is retained, so no Taylor truncation estimate is assumed. Adjacent-segment jumps and the initial MPS approximation are separately charged. Square roots are rounded upward by integer arithmetic. Arbitrary-size integers are used.

## 5. Observable intervals

Let n=||Phi0||^2, computed exactly from the supplied rational MPS. The exact evolution has the same norm. If ||Phi-q||<=eta and O is Hermitian, then

\[
\left|\frac{\Phi^\dagger O\Phi-q^\dagger Oq}{n}\right|
\le
\frac{\|O\|\,\eta(\sqrt n+\|q\|)}{n}.
\]

This follows by expanding Phi=q+error in the difference or using the two-sided product factorization. For D, ||D||<=2. The numerator and ||q||^2 come from exact small rational contractions. The approximation is **not** silently assumed normalized.

The population contrast refers to canonical orbital occupations, not measured real-space charge transport or a superconducting/battery property.

## 6. Uniform waveform, initialization, and noise coverage

### Waveform uncertainty

If an actual waveform satisfies

\[
A=\int_0^T(|u-u_0|+|v-v_0|)dt,
\]

then ||delta K||<=2(|delta u|+|delta v|). The normalized final state differs by at most 2A, and its D expectation by at most 8A. For a pointwise bound |delta u|,|delta v|<=a, this becomes 16Ta.

This covers every integrable waveform satisfying the bound, including nonsmooth changes. No assertion that D and W commute is needed. It is not a claim that every possible control is accurately represented by the small model.

### Initialization uncertainty

Let the actual normalized initial density matrix be within trace distance r0 of the specified normalized MPS, in the same balanced-spin sector. Unitary evolution, and the noisy evolution below, cannot increase that trace distance. Since ||D||<=2, the additional observable error is at most 4r0. The trace-distance condition is a mathematical assumption, not an experimentally measured preparation fidelity.

### Local phase-flip dephasing

Consider

\[
\dot\rho=-i[K(t),\rho]+\sum_j\gamma_j(t)(Z_j\rho Z_j-\rho),\quad
Z_j=I-2n_j,\quad\gamma_j\ge0.
\]

Every Zj is unitary and preserves the accepted sector. If

\[
\Gamma=\int_0^T\sum_j\gamma_j(t)dt,
\]

the trace distance to the corresponding noiseless state is at most Gamma. One proof uses trace-norm contraction and ||Zj rho Zj-rho||_1<=2. Equivalently, in a phase-flip jump representation the probability of any jump is 1-exp(-Gamma)<=Gamma. The D allowance is therefore at most 4Gamma.

This is a particular bounded Lindblad noise model, **not** arbitrary molecular thermalization, solvent interaction, or all environmental noise.

### Combined robust interval

If the nominal interval is [a,b], then

\[
\boxed{\langle D(T)\rangle\in
[a-16Ta_u-4r_0-4\Gamma,\;
b+16Ta_u+4r_0+4\Gamma].}
\]

The saved controls both certify D(T)<=-3/5 including their declared nonzero waveform, initialization, and dephasing budgets.

## 7. Numerical discovery and what it costs

All proposal work runs on ordinary CPU numerical libraries. It builds the full balanced-spin Hamiltonian (903,620 stored nonzeros), expands the initial MPS, samples numerical full-model trajectories, and computes POD/SVD bases. No new ground-state certificate is needed for a dynamical claim about a specified initial state.

An initial multi-control snapshot model failed on a different switching history even with 128 retained coordinates: its numerical integrated-defect estimate exceeded one. The retained successful bases are intervention-adapted, not validated as universal controllable subspaces. The short-pulse optimizer is local L-BFGS-B; the long pulse was selected from a finite exploratory scan. Neither is a proved global optimum or a new general optimal-control algorithm.

The final acceptance uses J coefficients on a 2^-32 grid, model matrices on 2^-40, and polynomial coefficients on 2^-48. The exact checks include the resulting effects instead of treating these grids as accuracy estimates.

A 24-coordinate short intervention and 64-coordinate longer intervention pass. The 8-coordinate version of the short intervention does not establish the requested target, although an independent numerical full calculation suggests the actual nominal pulse achieves it. Refusal means insufficient proof, not a proof of physical impossibility.

## 8. Limits and interpretation

The result adds model-level, finite-time, robust reachability and certified reduced dynamics. It does not provide a universal golden algorithm, cheap many-body discovery, experimental controls, exact ground-state control, synthesis reachability, a new physical law, room-temperature superconductivity, or a new H12 result. The electronic model has fixed nuclei and finite basis.

The main remaining computational challenge is to discover and validate equally useful embeddings directly in compact operator/tensor form rather than using the enumerated snapshot stage. This implementation does not establish that such an inexpensive embedding exists for a broad family.

The identities are explicit and computations independently integer checked, but they are not a formal theorem-prover proof and do not imply error-free software.
