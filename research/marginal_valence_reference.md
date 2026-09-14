# Valence reference resolves the finite H6 complement-gap blocker

Follow-up: [the charge-operator experiment](marginal_charge_operator_progress.md) replaces the expanded metric with a 24-factor one-site/two-site occupation recipe and retains both complete energy proofs. It also brackets the best possible bound of the simplest termwise occupation envelope, showing why that implicit route cannot meet the H6 gap threshold. Full determinant-row coverage remains explicit.

The localized H6 calculation now has a complete energy certificate with **20 retained determinants, 11 response directions, no dense complement factor, and width7.416980322924775e−12 Hartree**. A hopping perturbation of strength1/50 transfers using the same reference gap,14 newly discovered response directions, and width3.9293895949815965e−12 Ha. Both pass independent exact replay.

The decisive change is reference selection. The earlier localized P32 included only **5 of the20 singly occupied configurations**, together with27 ionic configurations. Its Q therefore retained most of the spin manifold, and its numerical Q floor was−6.311577514960223 Ha. Retaining all20 singly occupied configurations moves the numerical Q floor to−5.993874956660395 Ha. An exact weighted diagonal-dominance proof alone certifies Q≥−6.131369418195366 I for the original rational localized Hamiltonian. The previous reference target−6.26489910104 is now passed with no added three/four-coordinate atoms.

The old fixed-P cone obstructions remain valid. This construction changes P and therefore changes QHQ. It resolves the molecular gap requirement by selecting the whole valence manifold, rather than proving a stronger decomposition of the old complement.

## The reference representation

For each spatial site, let

\[
s_i=n_{i\alpha}+n_{i\beta}-2n_{i\alpha}n_{i\beta},\qquad
P=\prod_i s_i,\qquad Q=I-P.
\]

The commuting factor s_i is exactly1 on single occupation and0 on empty/double occupation. At half filling, every Q determinant has both a hole and a doublon. For six sites with Nα=Nβ=3, P has binom(6,3)=20 configurations. Direct site combinations generate these20 states; no Hamiltonian-sector diagonalization selects them.

The initial physical upper witness is a product of adjacent-site spin singlets,

\[
\bigotimes_{j=0}^{2} (|\uparrow_{2j}\downarrow_{2j+1}\rangle-
|\downarrow_{2j}\uparrow_{2j+1}\rangle).
\]

It has8 integer-amplitude determinants. Exact CAR evaluation checks that total spin raising annihilates it. Twenty-four physical Krylov steps then refine the witness using the actual Hamiltonian. This is a trial-state construction; the final ground-energy proof still checks the whole spin-zero sector, including its higher-total-spin representations. It does not assume the ground state is a singlet.

The projector has a short product description, but the implementation still enumerates its20 basis states and all380 ionic states. Its degree grows with the number of sites. Neither a constant-dimensional marginal representation nor efficient general scaling follows.

## Complete proof chain

For manageable rational arithmetic, the constructor rounds the exactly rotated Hamiltonian coefficients to denominator10^14 and restores exact SU(2) symmetry with the existing algebraic projection. It independently checks the resulting commutators. The exact coefficient-norm difference to the input Hamiltonian is **3.948152624701112e−12 Ha** and is recomputed in final replay.

A numerical comparison eigenvector proposes380 positive integer weights. Exact complete row checks certify a nearby-Hamiltonian Q lower of−6.131369418193787 Ha. The gap proof is an empty-atom certificate: only the weighted DD residual is needed. The numerical380-dimensional comparison eigensolve remains part of construction; eliminating dense Q factorization does not eliminate this explicit proposal matrix or the state enumeration.

The upper witness has Rayleigh value−6.333058626232914 Ha for the nearby Hamiltonian and exact variance3.521907636031548e−14. The excitation threshold is the short rational number−6.323058626233 Ha. Fresh directional Schur construction needs11 directions to prove at most one eigenvalue below that threshold. Exact inertia, the valid Rayleigh witness, and the strict complement bound establish the Temple premises. The final original-input interval is

\[
E_0\in[-6.333058626240384,\;-6.333058626232967]\ \mathrm{Ha},
\]

where endpoints shown are decimal approximations; the certificate stores exact rationals. The width is7.416980322924775e−12 Ha. Construction takes52.78 seconds after input preparation, including upper refinement and exact checks.

The energy certificate is for the input rational localized Hamiltonian. Its exact orbital-rotation proof from the preceding experiment establishes isospectrality with the saved canonical spin-symmetric reference. These are electronic finite-basis model energies, not continuum/basis-set error or thermal-property guarantees.

## Transfer under a physical hopping perturbation

The new transfer adds

\[
\Delta=\frac1{50}\sum_\sigma(a^\dagger_{0\sigma}a_{1\sigma}+a^\dagger_{1\sigma}a_{0\sigma})
\]

between the first two **localized spatial orbitals**. This differs from the earlier perturbation applied between canonical molecular orbitals. Two exact CAR-square identities bound the norm by1/25, attained by the hopping extremal eigenspaces, and the reference Q proof transfers with the full0.04-Ha scalar shift. The valence reference supplies enough gap to absorb that shift; it does not evade or contradict the rank/multiplicity theorem that prevented reducing the norm.

The physical reference witness seeds a new24-step upper refinement. All14 response directions are discovered on the perturbed Hamiltonian; no old response prefix is imported. Exact replay gives

\[
E_0(H+\Delta)\in[-6.305079525073250,\;-6.305079525069321]\ \mathrm{Ha},
\]

with exact interval width3.9293895949815965e−12 Ha. The Q component remains380-dimensional, and the largest dense complement factor dimension is0. Construction after source replay takes67.35 seconds. Final replay reconstructs the reference atom/DD proof, recomputes the norm identities and original-H error, and checks the actual perturbed response and witness.

## Charge-pattern weights

Independent numerical feasibility searches tie weights to progressively richer occupation features. All still inspect the explicit380 Q rows.

| Weight feature | Number of classes | Achieved numerical row bound, Ha |
|---|---:|---:|
| Doublon count | 3 | −6.522266368428464 |
| Count and absolute charge dipole | 19 | −6.338734231137632 |
| Charge pattern modulo reflection | 73 | −6.192147763127332 |
| Full charge pattern | 140 | −6.192224753447251 |

These are achieved proposals, not certified optima. Finite LP feasibility tolerances and the chosen near-boundary weights prevent interpreting the small ordering reversal between73 and140 classes as a cone-inclusion result.

The73-class proposal rounds to positive integers, and exact expansion into every physical Q row verifies an actual DD floor **−6.192147763123858 Ha**. The exported short threshold is **−6.2 Ha**. Reflection only ties the chosen weights; no reflection symmetry of the rational Hamiltonian is assumed. The same charge pattern receives the same weight regardless of spin labels.

The earlier11 response directions fail the exact inertia gate with this weaker gap, and that failure is retained explicitly. A fresh discovery needs **13 directions** and recovers the same7.416980322924775e−12-Ha energy interval. The73-weight reference also transfers under the same localized hopping:17 fresh response directions recover the same3.9293895949815965e−12-Ha perturbed interval, with no dense Q factor. Thus independent per-determinant weights are unnecessary for both finite proofs. `metric_rule.json` stores the73 generating weights; the present certificate format still includes their380 expanded values, so serialized proof compression and implicit row verification remain further work.

## Implementation and verification

`marginal_valence_states.py` generates the valence states and dimer singlet. `marginal_valence_reference.py` constructs the complete proof directly from the input Hamiltonian and sector labels. It imports no FCI vector or inherited reference/response basis.

`spin_gap` now accepts an explicitly selected atom/DD proof family both directly and inside a transferred reference. It rebuilds all blocks from the caller's H and P and retains the exact physical-state, coverage, residual, and positive-denominator gates. Conflicting direct proof families are rejected. The legacy factor path remains available. Existing transfer construction now accepts an independently verified Temple source and carries the atom/DD proof through the same norm-shift and directional-response machinery.

The full **302-test suite passes in285.536 seconds**. New tests cover exact single-occupation support, nonzero singlet spin raising, refusal at reference-budget overflow, input quantization budgets, a complete small valence energy interval, and actual atom-gap transfer. Independent `python -S` replay verifies the energy intervals, their physical upper witnesses, and the exact charge-gap certificates.

Artifacts are under `results/marginal_h6/valence_reference_seed`, `valence_temple`, `valence_transfer_1_50`, `valence_charge_classes`, and `valence_class_transfer_1_50`. The generated charge-class proposal and the rejected inherited response are labeled separately from accepted certificates.

## Remaining scope

The finite H6 complement-gap and dense-factor requirements are resolved for this valence reference, including the stated hopping transfer. The remaining scaling problem is concrete: certify the ionic sector from occupation rules and discover/compress the spin response without enumerating every configuration. Current replay still touches400 spin-sector determinants. At eight half-filled sites the valence manifold already contains70 spin-zero determinants, exceeding the current32-state reference cap; the helper explicitly refuses instead of truncating it silently. Class counts, response dimension, and rational bit growth remain uncontrolled as size increases.

The tested reference applies to this half-filled localized minimal-basis chain. General molecular reference discovery, marginal representability, metallic/thermal behavior, and material synthesis remain outside the proved results.
