# Certified many-body energy bounds

For pairwise anticommuting Hermitian Pauli strings, `(sum c_i P_i)^2 = sum(c_i^2) I`, so each clique gives the rigorous lower bound `-sqrt(sum(c_i^2))`. The implementation uses integer square roots and outward rational intervals; floating-point eigensolvers are not part of the proof.

This narrow certificate respects the complexity boundary in Liu, Christandl, and Verstraete, [N-representability is QMA-complete](https://arxiv.org/abs/quant-ph/0609125), which rules out silently claiming a universal efficient method. It is a small auditable counterpart to the convex reduced-state bounds of Kull, Schuch, Dive, and Navascués, [Lower Bounding Ground-State Energies of Local Hamiltonians Through the Renormalization Group](https://arxiv.org/abs/2212.03014).

The upper bound is an exact computational-basis expectation. Identity terms are exact. Unspecified terms are charged by the triangle inequality as the sum of absolute omitted coefficients, so an incomplete learned Hamiltonian is never presented as unconditional. These are total-energy bounds; energy-density claims require separate normalization and have different scaling.

The #9 interface is operator discovery: a later abstraction layer can propose reusable Pauli groups, while this checker rejects missing terms, invalid anticommutation, bad square-root intervals, or altered bounds. Greedy grouping is only a search heuristic. Exact diagonalization may validate small cases externally but is not a certificate.

The missing lemma for stronger prediction is a quantitative relation between operator/cluster size and certified energy-density error for a declared physical class. Failure cases include frustration, critical long-range correlations, and closures that break after adding an interaction or intervention.
