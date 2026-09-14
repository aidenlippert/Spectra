# Sources: interference routing and compressed lower-bound certificates

Scope: primary papers and exact matrix consequences relevant to the proposal. These sources support antecedents and limitations; they do not establish the proposed fermionic compression mechanism.

## Signed Laplacians and effective resistance

* **Chen, Wang, et al., “Characterizing the positive semidefiniteness of signed Laplacians via Effective Resistances,” IEEE CDC 2016.** [IEEE record](https://ieeexplore.ieee.org/document/7798396), [author PDF](https://eeqiu.people.ust.hk/wp-content/uploads/2021/09/Characterizing-the-Positive-Semidefiniteness-of-Signed-Laplacians-via-Effective-Resistances.pdf). For a connected signed graph, the signed Laplacian is PSD with a simple zero eigenvalue iff the generalized effective-resistance matrix associated with a spanning forest of negative edges is positive definite. Their Theorem 3 separates positive and negative subnetworks: in their notation, PSD/corank one requires `Gamma_F^+ < -Gamma_F^-` (with sign conventions defined in the paper). This directly supports a compressed *edge-space* test, but the matrix dimension is the number of independent negative edges, which can still be exponential in a configuration graph.

* **Zelazo and Buerger, “On the Definiteness of the Weighted Laplacian and its Connection to Effective Resistance,” IEEE CDC 2014.** [author PDF](https://connect-lab-technion.github.io/Publications/Zelazo_CDC14.pdf), Theorem III.3. With one negative edge `e-=(u,v)` and connected positive subgraph `G+`, the signed Laplacian is PSD iff `|w(e-)| <= 1/R_uv(G+)`, equivalently `w(e-) >= -1/R_uv(G+)`. Strict inequality gives corank one (equality is the threshold and may add a zero mode). For multiple negative edges, scalar tests require the paper's isolation hypotheses; shared cycles require retaining cross terms in a resistance matrix.

* **Song, Hill, Liu, “On Extension of Effective Resistance with Application to Graph Laplacian Definiteness and Power Network Stability,” 2019.** [arXiv](https://arxiv.org/abs/1906.07632). Extends resistance/conductance to disjoint node sets and gives PSD/corank-one criteria in terms of effective conductances. It also bounds the number of negative Laplacian eigenvalues by the number of negative edges. This is relevant to multiport compression, but does not provide a polynomial-size representation for an implicitly exponential graph.

## Routing inequality: exact algebra and what it does (and does not) prove

For an oriented incidence split `B_+ , B_-` and positive diagonal magnitudes `W_+ , W_-`, if a matrix `F` satisfies `B_+ F = B_-`, then

`B_- W_- B_-^T = B_+ F W_- F^T B_+^T`.

Consequently, the check `F W_- F^T <= W_+` is sufficient for `L_+ - L_- >= 0`. If `W_+` and `W_-` are positive diagonal (hence invertible), the equivalent inverse-form inequality is `F^T W_+^{-1} F <= W_-^{-1}`; `F` may be rectangular and need not be invertible. The equivalence follows by equality of the nonzero singular values of `W_+^{-1/2} F W_-^{1/2}` and its transpose, or by the standard contraction criterion. This is a direct congruence/Schur-complement calculation and is consistent with the multiport effective-resistance criteria above. It is not, by itself, a theorem that a small `F` exists for a fermionic configuration graph. The identity `B_+F=B_-` must hold globally, including orientations and cycle-space compatibility; checking sampled configurations cannot establish it.

## Grounding, diagonal slack, and mixed-sign diagonals

For an ordinary unsigned weighted graph, a grounded Laplacian is a proper principal submatrix of the PSD Laplacian and is positive definite when every connected component is grounded. See **Shao and Pan, “Distributed Stabilization of Signed Networks via Self-loop Compensation,” 2021**, [arXiv](https://arxiv.org/abs/2109.12555), which studies diagonal self-loop compensation for signed Laplacians and gives quantitative sufficient conditions for restoring stability/PSD.

The safe general statement for a proposed diagonal correction `S` is simply: certify `L_signed + S >= 0` directly, or certify a Schur complement after separating grounded variables. Nonnegative diagonal slack can rescue an otherwise indefinite signed Laplacian; arbitrary mixed-sign diagonal terms cannot be treated as harmless grounding. A negative diagonal entry alone does **not** imply that the total matrix has a negative Rayleigh quotient: edge terms may dominate it. A negative diagonal can force indefiniteness when a test vector makes the remaining quadratic-form terms vanish (for example, an isolated coordinate), but this requires checking the full quadratic form. Thus “grounding” results do not automatically cover `diag(phi_x^2(ell_x-b))` when some local-energy residuals are negative.

## SOS / relative-residual certificates

* **Li and Lu, “Quantum variational embedding for ground-state energy problems: sum of squares and cluster selection,” 2023.** [arXiv](https://arxiv.org/abs/2305.18571). Constructs a quantum SOS SDP hierarchy giving lower bounds on ground-state energy and relates it to RDM and embedding methods. This supports the established status of operator-SOS lower bounds and the use of structured cluster restrictions, while numerical tightness depends on hierarchy level and cluster choice.

* **Hastings, “Limitations and Separations in the Quantum Sum-of-squares, and the Quantum Knapsack Problem,” 2024.** [arXiv](https://arxiv.org/abs/2402.14752). Gives an explicit limitation: a restricted degree-4 Majorana SOS fragment can fail to give the correct order of magnitude for SYK ground energy. Therefore an approximate annihilator or a low-degree SOS fit is not a certificate unless all algebraic relations and the residual operator bound are verified.

For any exact decomposition `H-bI = S + Z_N + R`, with `S=sum A_a^dagger A_a` and `Z_N` zero on the target particle sector, the lower bound is `E_0 >= b - ||R||` by the variational principle. The relative alternative is `R >= -eta S-epsilon I` with `0<=eta<=1`, which gives `H-bI >= (1-eta)S-epsilon I` on that sector. This operator inequality must be certified throughout the sector; an expectation on a trial state does not imply a lower bound. Factored operator programs likewise require exact identities or rigorous enclosures throughout their declared domain.

## Implications for the proposal

1. The signed-Laplacian routing idea has established mathematical predecessors (multiport passivity/effective resistance), so it should be presented as a proposed application/compression of those criteria.
2. A small transfer matrix is plausible as a *certificate format* when a global incidence factorization and a PSD matrix inequality are supplied. No source found establishes that molecular configuration-space graphs admit polynomially many such templates.
3. Diagonal local-energy slack is an independent PSD burden. Mixed-sign residuals require a joint certificate; “grounding” intuition is valid only for the appropriate nonnegative slack/principal-submatrix setting.
4. SOS and relative-residual language is established, but restricted hierarchies have documented failure modes. Any claimed `0.0016`-Hartree gap therefore needs separately rigorous trial-state upper and operator lower bounds, with numerical roundoff controlled.

## CAR local rule and spin-chain antecedent

The displayed `K_pq` is a Jordan–Wigner-dressed fermionic hopping/projector operator. Its three-mode positivity calculation is an elementary local block check, and the coefficient threshold `t <= ab/(a+b)` is the harmonic-mean condition obtained from that block. The closest standard physics antecedent is the frustrated ferromagnetic spin-1/2 `J1-J2` Heisenberg chain: **Richter, Haertel, Ihle, and Drechsler, “Thermodynamics of the frustrated ferromagnetic spin-1/2 Heisenberg chain,” 2008/2009**, [arXiv](https://arxiv.org/abs/0811.3549). That paper studies ferromagnetic nearest-neighbor and antiferromagnetic next-nearest-neighbor couplings and identifies the `J2=|J1|/4` boundary for the ferromagnetic ground state. This is an antecedent for the coupling pattern and threshold scale, not evidence that the CAR construction is a new exactly solvable chemical Hamiltonian.

For the **spinless modes used here**, identify each occupation basis `{0,1}` with a two-state pseudospin. The explicitly included parity string cancels the Jordan–Wigner hopping sign. Directly in this tensor-product basis, `K_pq=I-SWAP_pq=2(1/4 I-S_p·S_q)` on the entire Fock space. No additional singly occupied physical-site projection is needed: these pseudospins represent mode occupation, not electron spin. Thus the uniform chain control is known frustrated ferromagnet physics, not a new solvable family. Fixed particle number selects fixed pseudospin magnetization. The open-boundary allocation and inhomogeneous couplings in the implementation still require their displayed local budget proof. Adding diagonal finite-range `V` requires an exact sector minimum and trial expectation; the implemented DP covers only the explicitly specified one-site fields and nearest-neighbor density potential.
