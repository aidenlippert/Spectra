# H6 compression: sharper response and exact spin reduction

Two finite costs have been reduced. Raising the certified complement threshold from -6.323058626212 to **-6.265 Ha** reduces the successful response basis from 32 to **28 directions**. Separately, an exact spin reduction lowers the final proof's explicitly represented configurations from 924 to **400**, while still bounding the original Hamiltonian's global fixed-number ground energy. The combined certificate has width **1.941992366627265e-10 Hartree**.

The reduction artifacts described here inherit a previously constructed full-Q factor-response certificate. The subsequent [direct spin constructor](marginal_h6_direct_spin.md) removes that dependency: from the Hamiltonian alone it selects P, builds its own upper, discovers and factors Q blocks, and certifies width1.93199e-10 Ha using400 configurations. Both sector dimensions still grow exponentially at fixed filling.

## Sharpening the complement threshold

The existing exact Q partition is reused, but every block is refactored for the new threshold and independently checked. All 12 blocks pass at integer scale 10^12. The smallest residual-certified positive margin is approximately 1.1248204479e-7 Ha. The new complement construction takes about 3.57 seconds on the recorded run.

The optimized response construction starts afresh with this stronger scalar threshold and the same independently evaluated Krylov upper. Its progress is:

| Response directions | Certified width (Ha) |
|---:|---:|
| 16 | 0.04604548909119924 |
| 24 | 0.0029483485591992368 |
| 25 | 0.00024084371419923667 |
| 26 | 5.490390199236662e-6 |
| 27 | 3.194561992366627e-7 |
| 28 | 1.401992366627265e-10 |

Construction took about 211.31 seconds, including repeated exact acceptance checks. Replay still covers all 924 determinants on this unreduced path. The earlier exact 31-direction lower bound was proved at a different complement threshold and does not constrain this sharper construction. The diagnostic negative count of 21 at the new threshold has not been exported as an exact inertia proof, and 28 is not asserted to be minimal.

## Why spin symmetry permits a global reduction

For paired alpha/beta spin orbitals, define the exact CAR operators

\[
S_+=\sum_p a^\dagger_{p\alpha}a_{p\beta},\qquad
S_-=S_+^\dagger,\qquad
S_z=\tfrac12\sum_p(n_{p\alpha}-n_{p\beta}).
\]

A number-conserving Hermitian Hamiltonian commuting with S_+ and S_z also commutes with S_-. Its energy eigenspaces are SU(2) representations. For even electron number, every spin irrep has integer spin and contains an S_z=0 state. Therefore

\[
E_0(H\mid N)=E_0(H\mid N,S_z=0).
\]

This includes high-spin ground states: it does not assume the ground state is a singlet. For H6, the right-hand sector has dimension binomial(6,3)^2=400 instead of binomial(12,6)=924.

The reduction applies to the global Hamiltonian before choosing P and Q. The retained P32 need not be spin adapted. Its complement is then taken within the 400-dimensional spin-zero sector. No assertion that the old full-sector Q projector commutes with spin is used.

## Handling coefficient rounding exactly

The original rational export does not commute exactly with S_+: its commutator has 108 nonzero CAR coefficients with coefficient L1 sum 1.08e-10. Approximate commutation alone is not accepted as a symmetry certificate.

The proposer constructs a nearby exactly spin-symmetric operator. On a number-conserving Hamiltonian with at most two-body terms, the adjoint spin representation has ranks at most 2. Define

\[
\mathcal C(X)=[S_z,[S_z,X]]+
\tfrac12\bigl([S_+,[S_-,X]]+[S_-,[S_+,X]]\bigr).
\]

Its rank-0/1/2 eigenvalues are 0,2,6, so the scalar projection is

\[
H_s=(I-\mathcal C/2)(I-\mathcal C/6)H.
\]

This polynomial is only a proposal mechanism. Replay independently verifies Hermiticity, number conservation, the even-particle hypothesis, and the exact commutators [H_s,S_+]=[H_s,S_z]=0. It does not trust a claimed projection or a small commutator residual.

Every CAR monomial has operator norm at most one. Re-expanding H-H_s therefore gives the rigorous full-space perturbation bound

\[
\|H-H_s\|\le\epsilon=
\sum_w|h_w-(h_s)_w|
=\frac{27}{500000000000}=5.4\times10^{-11}\ \mathrm{Ha}.
\]

If the spin-zero Schur proof certifies H_s>=L_s, the original global ground energy satisfies E_0(H)>=L_s-epsilon. The upper witness is separately evaluated on the original H, so no upward epsilon shift is needed for that endpoint.

## Exact reduced replay and its costs

The spin-zero oracle validates every state against both N and N_alpha=N_beta=N/2. Complete factor coverage now requires exactly 368 Q determinants. Only the two original blocks of dimensions 200 and 168 remain. Replay checks their complete coverage, absence of interblock Q coupling, and positive factor residual margins against H_s itself. Factor entries are rescaled without changing their rational values where the symmetrized coefficients require a larger common denominator.

The same sparse response basis and retained Schur check are reconstructed on H_s. The two reduced certificates give:

| Inherited response | Original-H certified width (Ha) | Certificate bytes |
|---:|---:|---:|
| 32 directions, original threshold | 1.761992366627265e-10 | 1,665,178 |
| 28 directions, sharper threshold | 1.941992366627265e-10 | 1,643,665 |

Both reduced replays use 400 distinct determinant sources and reference 400 determinants. More precisely, they evaluate 400 source actions under H_s and 200 under H for the independent upper, for 600 operator-specific source evaluations. These are not 400 evaluations of a single operator. The largest factor remains dimension 200.

The recorded reduction/export runs took about 11.39 and 10.21 seconds. They load an existing full-Q certificate, select its spin-zero blocks, rescale inherited factors, and perform exact replay. These timings exclude discovery of the inherited full-Q factors and response directions. No end-to-end reduced discovery speedup is claimed.

## Remaining work

The subsequent direct constructor now generates the reference, upper, complement partition, and response in the spin-zero problem from the bare Hamiltonian, without first generating a924-state proof. A spin-preserving perturbation can still merge the small complement blocks beyond the current factor budget. The harder requirement remains: control the growth of the 400-state sector, factor sizes, reference size, response dimension, and integer precision as the molecule grows.

Spin reduction supplies a sound structural improvement; it does not characterize the full cone of physical marginals or collapse generic many-body complexity. Spin-dependent Hamiltonians require the explicit perturbation transfer or a different symmetry argument. The present checker deliberately rejects odd-particle and unpaired-spin-orbital variants.

## Independent checks

```sh
python -S -m experiments.marginal_h6_complement --verify-factor results/marginal_h6/complement_sharp/certificate.json
python -S -m experiments.marginal_sparse_response --verify results/marginal_h6/complement_sharp_response/certificate.json
python -S -m experiments.marginal_spin_reduction --verify results/marginal_h6/spin_reduced/certificate.json
python -S -m experiments.marginal_spin_reduction --verify results/marginal_h6/spin_reduced_sharp/certificate.json
```

All four independent replays passed without site packages. Tests verify exact spin projection, rejection of approximate symmetry, a high-spin ferromagnetic ground state represented in S_z=0, a non-spin-adapted retained reference, complete reduced-sector coverage, the nonempty-complement gate, and a spin-breaking field whose true ground energy moves outside S_z=0. That last case checks that the recomputed perturbation bound, rather than an unchecked symmetry assumption, protects the global lower endpoint. False sharper thresholds and overwritten export paths are rejected.
