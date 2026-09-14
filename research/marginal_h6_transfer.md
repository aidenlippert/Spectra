# H6 transfer: exact reference obstruction and sparse upper adaptation

The larger molecular test separates two blockers that were hidden by the H4 examples. The selected H6 reference is provably inadequate for an accurate upper witness. A sparse physical Krylov refinement escapes that reference and sharply improves the upper. The complement row-bound certificate fails. The subsequent [sign-preserving complement and response proof](marginal_h6_complement.md) resolves that lower-bound blocker and certifies an H6 interval of width **1.221992366627265e-10 Ha**, while explicitly covering all924 determinants. This document records the upper adaptation and the obstruction that led to that proof.

## Fixture and independent checks

The fixture is a linear chain of six hydrogen atoms at 1.4 Angstrom spacing, in STO-3G with RHF canonical orbitals. It has12 spin orbitals,6 electrons, and a full fixed-number sector of924 determinants. PySCF2.14.0 supplies numerical integrals. The normal-ordered electronic Hamiltonian contains918 nonzero rational terms after Hermitian symmetrization and rounding to denominator10^12.

The exporter records the exact coefficient L1 difference from the canonical binary-float input, including discarded small coefficients and symmetrization. That quantity bounds the corresponding operator perturbation; it does not certify the numerical integral algorithm, basis error, or physical-model error. Nuclear repulsion is recorded separately from all electronic energies below.

Independent checks establish exact Hermiticity, total and separate alpha/beta number conservation, and agreement of the RHF determinant's CAR diagonal with the numerical RHF electronic energy. The independent FCI control uses the same spatial integrals in the alpha3/beta3 sector. Its CI vector is converted from alpha-then-beta ordering into interleaved spin-orbital ordering with the exact fermionic permutation sign, rounded into integers, and evaluated by exact CAR action.

The rational reference upper is approximately **-6.333058626233001 Ha**; the numerical FCI electronic value is approximately -6.333058626232415 Ha. Their difference is about5.86e-13 Ha. The FCI control is explicitly separate from the sparse constructor and is not a global lower certificate. Its integer witness has200 nonzero amplitudes.

## What the unchanged retained budget exposes

The sparse selector receives only the Hamiltonian and particle/mode labels. With its existing P32 cap, residual-first selection ends without an accepted complement-gap certificate. It evaluates37 source actions and references355 determinants during that failed attempt. These counts describe the selection attempt, not complete sector coverage.

Its retained upper is approximately -6.303696584750138 Ha. Comparing that number with a better upper alone does not measure its error exactly. We therefore exported a stronger statement using an exact retained-space floor.

Let ell be the checked rational floor and u_ref the separately evaluated variational upper. Exact LDL verifies

\[
H_{PP}-\ell I\succ0,\qquad \ell-u_{\rm ref}
=0.02936204138200049\ldots\ {
m Ha}>0.
\]

For every normalized psi entirely in P,

\[
\langle\psi,H\psi\rangle-E_0
\ge\ell-u_{\rm ref}>0.029362041382\ {
m Ha},
\]

because E0<=u_ref. This proves that **no upper witness confined to this selected P32 can support a1e-7-Hartree ground-energy interval**, regardless of how accurate the lower response becomes. It does not prove that every32-determinant reference is inadequate, or that a broader implicit reference needs the same number of coordinates.

The obstruction is independently replayable using only a32-dimensional retained positivity check and the exact reference Rayleigh quotient. Its numerical FCI origin is unnecessary for the validity of the upper once the integer vector is supplied.

## Escape the retained space through physical Hamiltonian actions

The new sparse upper proposer starts from the failed selector's own integer witness. It builds a small sequence of physical Hamiltonian-action directions, with numerical orthogonalization and integer rounding. Projected generalized eigensolves propose a new physical vector, which is rounded and checked by exact Rayleigh evaluation. The independent FCI witness is not supplied to this refinement.

Rounding means the proposal basis is not asserted to be an exact mathematical Krylov sequence. The accepted object is its final integer physical vector. Replay also computes its eigen-residual norm squared exactly, without interpreting that residual as a ground-state error bound.

| Upper construction | Support size | Difference from independent rational reference upper (Ha) | Exact squared eigen-residual, shown numerically (Ha²) |
|---|---:|---:|---:|
| Initial selected P32 | 32 or fewer | 0.029362041482862988 | Not used as an error certificate |
| 12 physical Krylov directions | 200 | 1.2285589362533887e-6 | 9.958749746084823e-7 |
| 20 physical Krylov directions | 200 | 2.1199724738023203e-11 | 1.9390701959838745e-11 |

The20-direction upper is approximately **-6.333058626211801 Ha**. Both refinements evaluate and reference200 determinants, rather than constructing the924-state sector matrix. Initial recorded construction times are about1.81 and3.96 seconds. These counts apply to the upper-refinement phase; they do not replace or hide the separate initial selector and FCI-control costs.

Agreement between two upper witnesses is not a lower bound on E0. Likewise a small eigen-residual locates proximity to some spectral value but does not identify the ground state without further information. Neither observation establishes the requested H6 ground-energy accuracy.

## The initial lower gate fails

Using the improved upper u20, the requested complement threshold is gamma=u20+0.01 Ha. Complete occupation-tree construction stops at a genuine complementary row, state bitmask3040, with

\[
g_Q(3040)=-6.420084816082\ {
m Ha},
\qquad \gamma-g_Q(3040)\approx0.09702618987019923\ {
m Ha}.
\]

This rejects the current sufficient Gershgorin row proof. It does **not** prove that lambda_min(H_QQ)<gamma, that the physical complement has no gap, or that no32-dimensional lower reference can work. The stored diagnostic includes the exact diagonal and complementary radius so that this approximation loss can be investigated independently.

This identified the need for a stronger complementary-energy proof or an adapted lower reference, while retaining the improved sparse upper. The subsequent factor proof supplies that missing complement bound. Blindly increasing response precision cannot bypass an unproved complement bound. Blindly reusing the old P-supported upper is now ruled out by the exact obstruction.

## Validation and scope

The fixture is generated by `.venv-molecule/bin/python -m experiments.marginal_h6_fixture`. The numerical generator explicitly requires converged RHF and FCI controls. The sparse construction, obstruction checker, and upper replays do not require PySCF.

Independent `python -S` replays verify the reference upper, the retained-space obstruction, and both new sparse uppers. Tests check all alpha/beta interleaving signs on a smaller complete basis, molecular number/spin conservation, exact RHF and FCI-derived upper evaluations, obstruction tampering, bare-Hamiltonian reproduction of the failed larger selector, the post-refinement row failure, exact residual arithmetic, and physical upper construction on an independently solvable example.

These results extend the investigation from70 to924 configurations and change the next action with a proved obstruction. These upper-only and obstruction artifacts do not themselves supply a global lower endpoint. The subsequent complement-response certificate does; a universal marginal-cone representation and an asymptotic efficiency theorem remain open. The exact200-state support on this symmetric chain is not evidence that generic perturbed molecules retain such a frontier.
