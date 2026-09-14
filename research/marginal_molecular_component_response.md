# Recovering the molecular complement response

The H4 square's scalar-Schur precision gap is now isolated and closed on the unperturbed fixture. Both existing 32-determinant retained spaces yield exact energy intervals of width **1.0025118312170343e-10 Hartree** when the full complement response is retained. Exact polynomial moments reproduce the same rational endpoints and Schur pivots as direct inverse solves.

This remains a finite construction on an explicitly assembled 70-state matrix. It establishes a response representation for these molecular tests, not a scalable marginal-cone characterization.

## What the controlled comparison resolves

The same Hamiltonian, retained space, and integer upper witness are used throughout. Only the complement approximation changes.

| Complement treatment | Coupling-selected P32 width (Ha) | Witness-selected P32 width (Ha) |
|---|---:|---:|
| One global Gershgorin scalar | 0.1304555622232512 | 0.28635766589125117 |
| One Gershgorin scalar per Q graph component | 0.10034571263025119 | 0.13599904928725118 |
| Exact-checked spectral lower per Q component, still scalar response | 0.10034571128625118 | 0.059136770454251184 |
| Exact inverse response per Q component | 1.0025118312170343e-10 | 1.0025118312170343e-10 |
| Exact moment response per Q component | 1.0025118312170343e-10 | 1.0025118312170343e-10 |

Better scalar energy thresholds alone do not recover the missing response. The exact inverse and moment variants preserve its action on the retained coupling directions. For both selections their complete rational Schur pivot lists agree exactly.

## Component bound and mandatory coverage

For exact connected components C_j of the induced Q graph, write D_j=H_CjCj and B_j=H_PCj. The Gershgorin bound gamma_j gives D_j>=gamma_j I. If every gamma_j>b, the sufficient lower certificate is

\[
H_{PP}-bI-\sum_j\frac{B_jB_j^T}{\gamma_j-b}\succ0.
\]

This dominates the old global scalar bound at fixed P, because each gamma_j is at least their global minimum and every B_j B_j^T is positive semidefinite. A Q-only component with B_j=0 still requires its gap check: it could contain a lower state even though it contributes no Schur correction.

The spectral variant proposes stronger gamma_j from small numerical eigensolves, then checks D_j-gamma_j I by exact LDL before using them. If the proposed value equals the original Gershgorin threshold, its analytic Gershgorin proof suffices. Numerical eigenvalues are not accepted as proof.

## Exact response from a checked polynomial

For a component coupling W_j=B_j^T, the generic exact recurrence routine discovers a monic polynomial f_j such that

\[
f_j(D_j)W_j=0.
\]

It checks the relation by substitution into the original operator powers. With q_b(x)=(f_j(x)-f_j(b))/(x-b),

\[
W_j^T(D_j-bI)^{-1}W_j
=-\frac{1}{f_j(b)}\sum_k q_{b,k}\,W_j^T D_j^k W_j.
\]

This positive self-energy is **subtracted** from H_PP-bI. Replay requires f_j(b) nonzero and all component gap checks, including uncoupled components. Only retained columns that actually couple to the component enter its moment matrices.

The coupling-selected square has Q component sizes 6,6,7,6,1,7,4,1 and response degrees **6,6,5,5,0,5,4,0**. The witness-selected space has sizes 6,12,6,1,6,6,1 and response degrees 6,12,6,1,5,5,1. Degree zero denotes no coupling, not permission to omit the component's energy bound.

The existing matched-reference recurrence was factored into a generic sparse operator-action routine and reused here. Moment replay does not invoke the direct inverse solver; a test patches that solver to throw and still replays the molecular moment certificate successfully. Both moment exports are about 46.9 KB, mostly the shared Hamiltonian and upper witness.

## Stressing exact component structure

Add the real, number-conserving hopping perturbation

\[
\varepsilon\sum_{i=0}^{6}(a_i^\dagger a_{i+1}+a_{i+1}^\dagger a_i).
\]

At epsilon=1/1000, the fixed P32 complement becomes one connected block of dimension 38. Its global Gershgorin lower still exceeds the old valid upper by 0.05846490449 Hartree. Nevertheless, direct recurrence discovery on the merged block exhausts the degree-12 budget. That is an observed compression failure, not an instability or absence of a physical state.

To retain the small reference blocks, partition Q using the old components and write the actual perturbed complement as H_QQ=C+R, where C contains its actual within-block entries. The exact bound

\[
\eta=\max_i\sum_j|R_{ij}|\ge\|R\|_2
\]

implies H_QQ-bI>=C-(b+eta)I. If every gamma_j>b+eta, inverse order gives the sufficient certificate

\[
H_{PP}-bI-\sum_j B_j\bigl(C_j-(b+\eta)I\bigr)^{-1}B_j^T\succ0.
\]

Moments use the actual perturbed within-block matrices and actual couplings. The response is exact for this block-diagonal surrogate; the interblock effect is bounded by eta. Replay verifies a complete, disjoint Q partition and recomputes eta rather than trusting a serialized norm claim.

| Hopping strength | Interblock norm bound | Block-moment width (Ha) | Full connected-inverse control width (Ha) |
|---|---:|---:|---:|
| 1/10000 | 7/10000 | 1.027610466794303e-8 | 1.0010466794302913e-10 |
| 1/1000 | 7/1000 | 1.2116563809476015e-5 | 1.0080947601473919e-10 |

Both stress cases retain P32 and use response degrees 6,6,7,6,1,7,4,0. Fresh integer upper witnesses are proposed by full 70-state numerical diagonalization and checked by exact Rayleigh evaluation. That upstream enumeration remains part of the method's cost.

The full-inverse controls isolate the wider intervals as a loss from the interblock norm approximation. The norm-shift approximation alone leaves the stronger perturbation above the 1e-7-Hartree target. The subsequent [directional interblock response](marginal_interblock_resolvent.md) closes this gap: second order gives 2.57809e-10 Ha with the same partition, and automatically selected blocks give 1.71809e-10 Ha. Both weaker-perturbation variants give 1.00105e-10 Ha.

## Verification and limits

All eight unperturbed component-method exports and all four connected-stress/control exports have independent `python -S` replays. The tests cover strict component gap checks, uncoupled low states, exact graph edges, false spectral thresholds, inverse/moment equality, recurrence exhaustion, changed Hamiltonians, partition omission/duplication, and the necessary interblock shift.

No bound here controls component growth, recurrence degree, rational bit complexity, or reference discovery for general molecules. The stress test explicitly demonstrates the fragility of exact graph components and the finite range of the norm-shift repair. Molecular errors outside the exported rational finite-basis electronic Hamiltonians remain separate.

The original component milestone passed **161 tests in 51.442 seconds**; current validation is recorded in `results/marginal_final_validation.json`. The exact-only replay commands use `python -S`; numerical proposal tests use the normal Python environment with NumPy.
