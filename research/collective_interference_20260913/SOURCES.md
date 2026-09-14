# Primary-source audit: grouped baths, bright modes, and collective bounds

## Bravyi–Gosset impurity precedent

**Bravyi and Gosset, “Complexity of quantum impurity problems,” CMP 356 (2017)**, [arXiv](https://arxiv.org/abs/1609.00735), studies a finite interacting impurity coupled to an otherwise quadratic free-fermion bath. They give a classical ground-energy approximation with runtime n³ exp[O(b³)] for additive error 2⁻ᵇ, and represent low-energy states by a controlled number of Gaussian states. The theorem concerns a quadratic bath plus an arbitrary Hamiltonian on O(1) impurity modes; it does not establish the band-endpoint certificate implemented in this pass.

## Bright/dark compression and what is exact

For a quadratic bath with one-particle matrix E and coupling matrix V, an SVD of V identifies rank(V) coupling directions (“bright” modes); their orthogonal complement is dark only when the bath one-particle Hamiltonian preserves that decomposition. Degenerate bath bands make this exact: within each degenerate eigenspace, rotate the basis so only the coupling-column span is bright and the remainder is dark. This is finite-dimensional linear algebra, not a new impurity theorem. For nondegenerate or merely grouped energies, the rotation generally creates bright–dark mixing; discarding dark modes then needs an inequality or another error bound.

Write the bath relative to a filled negative-energy reference as `H=E_filled+H_active+H_coupling+sum_{j in +}|e_j| n_j+sum_{j in -}|e_j|(1-n_j)`. Since each `n_j` and `1-n_j` is PSD, replacing each excitation coefficient by a lower or upper band endpoint gives a termwise operator order. With identical active Hamiltonian, couplings, and interactions on both sides, this yields `H_- <= H <= H_+` on the same Fock space; interactions cannot invalidate adding the same operator to both sides. The filled-reference convention is algebraic bookkeeping, not a ground-state occupation assumption. Only afterward should each degenerate bath eigenspace be rotated into bright and dark directions. No primary source located here gives this exact band-extrema construction, so it is recorded as the proposal's derivation rather than attributed prior art.

## Collective dark spectators and charging terms

If the interactions have the form F(N_active,N_bath), substitute N_bath = N_total − N_active on the fixed-number sector. N_active remains an operator: active–bright hybridization must be retained. Within each allowed total number of retained active-plus-bright particles, minimize the dark occupations by their degenerate energies and multiplicities. This is the exact elimination proved in PROOF.md. Additional dark-mode dispersion, noncollective interactions, or orbital-dependent couplings require their own argument; the envelope construction handles the specified dispersion before elimination.

The one-leg tensor-rank pitfall is narrow: a rank obstruction concerns the full-Fock interaction tensor's orbital support under orbital rotations. It does not account for fixed-number identities, rule out an operator sandwich, or preclude collectively minimized interacting spectators. The compact control itself has full one-leg rank before its number identity is used.

## Novelty boundary and limitations

Impurity reductions and bright-channel rotations have established antecedents. This pass supplies a particular proved and replayable band-extrema sandwich, exact filled-reference accounting, and collective fixed-number spectator minimization. The repository already contains implicit Dicke–Schur certificates. No first-in-literature claim is made, and a primary-source search is not a novelty proof. The Bravyi–Gosset theorem does not establish that generic molecular interactions fit this restricted collective grammar.
