# Adaptive quartic cut search

The bounded repaired run was completed on the matched ten-mode, five-particle model. It starts from the complete cubic SOS blocks and, at each step, adds eight rank-one quartic directions selected from the most negative eigenvalues of the current quartic moment matrices. The repaired run is stored in [`results/marginal_cut_search/cache_fixed`](../results/marginal_cut_search/cache_fixed).

All seven exported certificates replay through the exact interval checker. The physical upper endpoint is 3.28100669557401. The lower endpoints for iterations 0 through 6 are respectively:

```text
3.25715646060795  3.25715424337559  3.25715202576607  3.25715451304531
3.25715249805976  3.25714562369171  3.25714331790570
```

Thus the best width is the initial 0.02385023496606; the final width is 0.02386337766831. The cuts do not improve the energy bound. The selected moment violations also do not provide a useful progress metric: the worst eigenvalue moves from −1.031 to −0.427, then grows to −1.851. The exact certificates remain valid because their residual budgets are charged by the checker; solver status is `optimal_inaccurate` from iteration 2 onward.

The result diagnoses a selection failure, rather than a quartic impossibility. The violated directions are unconstrained higher moments chosen by spectral size. They need not overlap strongly with the Hamiltonian's energy gradient after the lower-degree affine and number-ideal constraints are imposed. Adding them can spend degrees of freedom on moment completion while leaving the energy-supporting face unchanged. A useful next selector must score a violated direction by its attainable energy effect (for example, a primal-dual sensitivity or a small reoptimization with the candidate cut), not by negative eigenvalue alone.

This run also confirms that the earlier interruption was an implementation cache bug: Boolean creation labels were colliding with integer labels in Python's word-product cache. The core canonicalization and adjoint paths now normalize those labels safely.

## Later controlled result

The seven-step failure above must not be generalized to all negative-eigenvector selection. A later symmetry-reduced run with 117 single cuts reached a certified width of 8.157088389144001e-7. Energy lookahead showed no advantage under that common cut budget. See `research/marginal_energy_selector_results.md`. The changed formulation and larger budget prevent attributing the improvement to only one change.
