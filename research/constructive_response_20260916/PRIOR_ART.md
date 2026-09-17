# Prior-art audit: constructive certified response compression

Date: 2026-09-16

## Scope

The proposed target is a direct procedure for a partition (P+Q=I) that constructs a compact operator enclosure for the Schur/Feshbach response

\[
\Sigma(E)=PHQ(QHQ-E)^{-1}QHP,
\]

with a certified Loewner error over an energy window, without first enumerating or solving the full (Q) many-body space, and with bounds that remain usable after repeated eliminations. I checked primary papers and abstracts/author pages linked below. I did not perform a comprehensive patent search, reproduce the algorithms, or claim that the proposed construction is new.

## Strongest precedents

### 1. Feshbach--Schur effective Hamiltonians

Dusson, Sigal, and Stamm give a self-contained Feshbach--Schur map, fixed-point formulation, and explicit eigenvalue/eigenvector estimates: [arXiv:2105.02058](https://arxiv.org/abs/2105.02058). The response is exactly the familiar Schur complement, and the reduced eigenproblem is nonlinear in (E). This establishes the reduction identity and rigorous spectral estimates. It does not provide a generally cheap, many-body, matrix-free construction of a two-sided Loewner enclosure for the response. The paper's applications are perturbative/discrete-spectrum examples, including helium-type ions.

### 2. Block Lanczos and resolvent quadrature

Block Lanczos has long been used to project resolvents onto a small Krylov block and represent them by matrix continued fractions: [The block Lanczos algorithm and the calculation of matrix resolvents](https://doi.org/10.1016/0010-4655(89)90155-0). Fenu, Martin, Reichel, and Rodriguez give block Gauss/anti-Gauss bounds for matrix-valued functions under stated conditions: [SIAM J. Matrix Anal. Appl. 34 (2013), 1655](https://doi.org/10.1137/120886261). More directly, Zimmerling, Druskin, and Simoncini's 2025 paper defines the block Gauss--Radau construction for (B^T\phi(A)B), proves matrix inequalities F_m\preceq F\preceq \widetilde F_m, monotonicity, and computable bounds for resolvents/Stieltjes functions: [Journal of Scientific Computing 103, 5](https://doi.org/10.1007/s10915-025-02799-z). These are genuinely Loewner-order results for the matrix-valued quantity under their symmetric positive-definite and spectral hypotheses.

This is the closest algorithmic precedent. A natural response approximation is to choose (B=QHP), apply block Lanczos to (QHQ), and evaluate (B^T(QHQ-E)^{-1}B) with Gauss/Gauss--Radau enclosures. The unresolved issue is the many-body one: obtaining the required Lanczos moments, a certified spectral interval, and a compact representation of the block-valued remainder without constructing the (Q) action or a basis whose cost scales like the full determinant space. Block quadrature bounds do not by themselves solve that construction problem.

### 3. Schrieffer--Wolff and local Schrieffer--Wolff

Bravyi, DiVincenzo, and Loss develop exact and perturbative Schrieffer--Wolff transformations, linked-cluster structure, and ground-energy precision bounds for local spin systems: [arXiv:1105.0675](https://arxiv.org/abs/1105.0675). The local-SW line gives locality-aware truncation estimates in perturbative regimes, while also identifying limitations of a local generator outside that regime. This is a strong precedent for eliminating high-energy sectors and retaining local effective interactions. It is not a nonperturbative, arbitrary-fermion, two-sided certified response enclosure whose construction remains compact through recursive elimination.

### 4. Density-matrix embedding theory

Knizia and Chan construct a small quantum bath for strongly coupled fragments and apply it to hydrogen rings/grids: [Density matrix embedding: A strong-coupling quantum embedding theory](https://doi.org/10.1021/ct301044e). Wouters et al. describe the practical bath construction and self-consistency in quantum chemistry: [A practical guide to DMET in quantum chemistry](https://doi.org/10.1021/acs.jctc.6b00316). DMET is directly relevant because it compresses environmental entanglement into a bath and can treat covalent coupling. However, its bath is generated from an approximate/reference state and its self-consistency is not an exact upper/lower enclosure of the eliminated Hamiltonian response. It therefore supplies a construction strategy and benchmark family, not the requested rigorous certificate.

### 5. Dynamical mean-field theory and impurity mappings

DMFT maps a lattice problem to a self-consistent impurity problem and is exact in its large-coordination/infinite-dimensional limit; finite-dimensional applications approximate nonlocal correlations. The primary review is Georges, Kotliar, Krauth, and Rozenberg, *Dynamical mean-field theory of strongly correlated fermion systems and the limit of infinite dimensions*, [Rev. Mod. Phys. 68, 13 (1996)](https://doi.org/10.1103/RevModPhys.68.13). DMFT demonstrates that a compact environment response can be powerful, but its exactness condition is a special limit and it does not furnish finite-system rigorous Loewner bounds for the discarded spatial correlations.

### 6. A provable impurity special case

Bravyi and Gosset prove a quasi-polynomial classical algorithm for a free-fermion bath coupled to an (O(1))-mode interacting impurity, with an explicit additive-energy guarantee and a Gaussian-state superposition representation: [Complexity of quantum impurity problems](https://arxiv.org/abs/1609.00735). This is a valuable positive precedent for compact environment descriptions, but its quadratic bath and constant-size impurity are precisely the restrictions that the proposed general correlated-fermion response would need to exceed.

### 7. The overlap obstruction

The many-body overlap needed by a global low-rank or guiding-state compression can decay with system size through Anderson's orthogonality catastrophe. The original primary result is Anderson, [Infrared catastrophe in Fermi gases with local scattering potentials](https://doi.org/10.1103/PhysRevLett.18.1049); a later interacting treatment is Yamada and Yosida, [Orthogonality Catastrophe Due to Local Electron Interaction](https://doi.org/10.1143/PTP.59.1061). A modern embedding discussion explicitly connects local approximation errors to global near-orthogonality: [High Ground State Overlap via Quantum Embedding Methods](https://doi.org/10.1103/PRXLife.3.013003). This does not rule out response certification, because energy/observable bounds need not require global state overlap; it does rule out treating a bare global low-charge or mean-field projector as a generally stable compressed representation without an additional decay or response theorem.

### 8. Local Markov and recovery structure at finite temperature

Chen and Rouzé prove that Gibbs states of bounded-degree quantum Hamiltonians are locally Markovian at arbitrary temperature, with quasi-local recovery maps and exponential conditional-mutual-information decay for shielded regions: [Quantum Gibbs states are locally Markovian](https://arxiv.org/abs/2504.02208). Earlier work by Kato and Brandão proves related approximate-Markov/thermal equivalences in one dimension: [Quantum Approximate Markov Chains are Thermal](https://arxiv.org/abs/1609.06636). These results are relevant to response compression, but their outputs concern state recovery/local marginals and include dimensional, clustering, or preparation assumptions. They do not automatically produce a two-sided ground-energy Schur enclosure for a strongly correlated 2D fermion Hamiltonian.

### 9. Local patch lower bounds and certified equilibrium relaxations

The basic cluster lower bound for a decomposition (H=\sum_j h_j), with correct counting/normalization, is standard in tensor-network practice and certifies energy density while usually losing boundary and correlation information. A systematic primary result is [Certified algorithms for equilibrium states of local quantum Hamiltonians](https://doi.org/10.1038/s41467-024-51592-3), which formulates convex matrix-inequality relaxations for local thermal and ground-state observables. Thus “small local certificates” or a planar local-SDP energy bound alone would not be a breakthrough. The new claim would have to retain inter-cluster response with certified size-controlled error and compose it across reductions.

## Why scalar residual identities are insufficient

A scalar residual estimate for one vector (x), such as

\[
|x^T(\Sigma-\widehat\Sigma)x|\leq r,
\]

certifies only one quadratic form. It does not imply a Loewner enclosure on the whole boundary space. Even sampled vectors can miss a bad direction without a certified covering or operator-norm bound. To obtain ( -R\preceq\Sigma-\widehat\Sigma\preceq R), one needs a matrix-valued construction, a certified norm bound, or a complete positivity/dual argument. Block Gauss--Radau already supplies this matrix-valued Loewner mechanism once block Krylov data and a spectral interval are available. Finite global Schur scalar identities are therefore useful moment checks, but cannot by themselves improve on known block-quadrature certification or establish compression.

## What appears genuinely open

## Audit of the new residual theorem

The proposed block result is best described as a specialized a-posteriori lower-bound lemma for a self-adjoint block operator. Its ingredients are established: block Schur complements/Feshbach maps, completion of squares by triangular congruence, operator Riccati graph transforms, and residual-based eigenvalue estimates. Relevant primary operator-theory precedents include [Solvability of the operator Riccati equation in the Feshbach case](https://arxiv.org/abs/1712.05770) and [Reducing graph subspaces and strong solutions to operator Riccati equations](https://arxiv.org/abs/1307.6439). Kato--Temple/Weinstein-type a-posteriori bounds likewise use a trial Rayleigh value, a residual, and spectral separation; a primary modern treatment is [On Temple--Kato like inequalities and applications](https://arxiv.org/abs/math/0511408).

The specific bound

\[
H-eI=\begin{bmatrix}A&B^*\\B&D\end{bmatrix},\quad D\succeq\delta I,\quad
K=A-B^*X-X^*B+X^*DX\succeq0,
\]

with (R=B-DX), (R^*R\preceq\rho^2(I+X^*X)), and
(2\rho<\delta), implies

\[
H\succeq\left(e-\frac{\rho^2}{\delta-2\rho}\right)I
\]

by a direct completion-of-squares/congruence argument. I found no source stating this exact normalization and constant in this exact notation, but that is not evidence of priority: it may be a useful repackaging of standard Riccati residual estimates. Its real novelty, if any, must come from a new many-body construction of (X,K,R,\delta,\rho) without enumerating the eliminated sector, together with a scope theorem and composability—not from the inequality alone.

The theorem also has a clear limitation. It is a global operator bound on the boundary residual in the metric (I+X^*X). If (R^*R\preceq\rho^2(I+X^*X)) is established by forming a full boundary Gram matrix, the global many-body cost has merely moved into the residual check. If it is replaced by scalar tests or a few sampled vectors, the Loewner conclusion is lost. The useful experimental question is therefore whether the residual metric can itself be generated and certified from compact local contractions.

The potentially new statement is the conjunction, not any individual ingredient:

1. Given a finite-range or chemically structured interacting fermion Hamiltonian and a declared (P/Q) partition, construct the response enclosure directly from local operator data and a certified gap/interval.
2. Bound the omitted response in Loewner order, or in an observable-specific order strong enough to certify the target energy, with no hidden full-space enumeration.
3. Preserve a compact, composable representation after the induced interactions are fed into a second elimination.
4. Give an a priori or cheaply checkable tractability criterion based on physical input, rather than discovering after a full solve that the representation happened to be small.

I found no primary source above that establishes all four properties for generic strongly correlated molecular or finite-dimensional fermion Hamiltonians. This is a research gap, not evidence that the goal is possible. Worst-case hardness and known failures of generic tensor-network/Monte-Carlo methods mean a theorem must state a restricted class or an adaptive certificate/failure output. The first credible breakthrough would therefore be a theorem of the form: for a clearly defined class (for example, a verifiable local gap/decay condition), the constructor produces a response enclosure of size and error bounded independently or mildly in system size, and recursive composition preserves those bounds.

## Feasibility judgment

The most promising near-term route is a hybrid of block-Lanczos resolvent enclosures and local/overlapping operator blocks. Lanczos/Gauss--Radau supplies a mathematically exact way to turn a spectral enclosure into two-sided response bounds; locality or tensor structure must supply a way to generate the needed block moments without full determinant enumeration. The decisive experiment is therefore not another H-size point. It is whether the first few response moments can be generated and verified from compact local contractions on a fully coupled case, with a residual enclosure that remains useful after a second elimination.

Do not call a post-hoc low-rank fit, a DMET bath, a perturbative SW truncation, or a small Krylov projection a world-level breakthrough by itself. Each is established territory unless it comes with the missing nonperturbative, directly constructible, composable certificate and a scope theorem.

## Source limitations

I inspected the primary source records and available abstracts/pages for the works linked above. I did not inspect every cited reference inside them, run their codes, search patents, or verify claims against unpublished material. The novelty conclusion is consequently a scoped literature audit, not a priority search or proof of originality.
