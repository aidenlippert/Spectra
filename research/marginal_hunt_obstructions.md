# Adversarial audit: the physical-marginal sandwich

This memo attacks the proposed compiler at the points where a valid variational sandwich can silently be mistaken for an efficient one. The architecture remains useful, but its breakthrough claim must be stated as a falsifiable complexity and certification hypothesis.

## 1. The cone does not erase the hard part

For fixed particle number, the set of normalized 2-RDMs is a compact convex **body**; its conic hull is the corresponding cone. Minimizing a two-body Hamiltonian over the exact body is equivalent to the ground-state problem. Thus “only one direction matters” is an instance-wise observation, not a complexity escape. The fermionic (N)-representability decision problem has QMA-completeness results under standard promise and encoding assumptions (Liu, Christandl, Verstraete, 2007: https://arxiv.org/abs/quant-ph/0609125); this does not classify every chemically structured directional instance.

**Repair.** Define the target as a restricted physical family: basis size (M), particle number (N), locality/decay class, accuracy ε, and certificate size/degree bounded by explicit functions. The claim is then: for a stated chemistry/materials distribution, the directional support function admits short *verifiable* certificates with high coverage. Test worst-case adversarial instances separately.

## 2. Separation exists, but short separation need not

Hahn–Banach supplies a separating functional, with no bound on its description length, algebraic degree, sparsity, or numerical conditioning. A separating hyperplane for a 2-RDM is simply a 2-body operator inequality valid on the (N)-fermion sector. If all such inequalities had uniformly compact efficiently checkable descriptions, one would obtain an unexpectedly strong resolution of the known QMA-hard representability problem. The machine may discover a useful cut, but it cannot infer universal validity from the violated sample.

**Repair.** Separate discovery from proof. The oracle proposes an operator; a deterministic verifier checks a certificate in a predefined proof system. Record failures as a certificate-complexity curve: minimal support, SOS degree, coefficient bit length, and verifier runtime versus (M,N,ε). “Sparse” must be measured in a fixed orbital/operator basis and compared with dense baselines.

## 3. SOS rank is not SOS degree

The displayed identity allows cubic B_α containing three ladder operators. Its square generally contains six ladder operators (three-body terms), not merely two-body terms. CAR normal ordering does not generically cancel those terms. A small number (m) of squares therefore says little: low rank can coexist with high degree, and increasing degree enlarges the dictionary rapidly. At fixed degree the full dictionary has polynomial size in the number of modes. This is exactly why the (P,Q,G,T_1,T_2,…) hierarchy is indexed by operator degree/particle rank, not just matrix rank (Mazziotti, 2012: https://doi.org/10.1063/1.3681138).

**Repair.** Require an explicit projection Π₂ onto the 0-, 1-, and 2-body quotient modulo the number ideal. Projection alone is not an operator identity: the discarded ideal terms must be shown Hermitian (or bounded on the fixed-(N) sector) and vanish there. Prove

\[
  \Pi_2\!\left(\sum_\alpha B_\alpha^\dagger B_\alpha\right)=\widehat A-bI
\]

with every 3+-body coefficient certified zero before quotienting. Track both SOS degree (d) and rank (r); report “(r,d)-certificate,” never rank alone. A useful first experiment is to enumerate exact low-degree certificates for small Hubbard molecules and plot the Pareto frontier ((d,r,ε)).

### A bounded CAR sanity check

There is a clean exact baseline for implementation. For (B=\sum_{i<j}c_{ij}a_j a_i), the matrix

\[
  D_{ij,kl}=\langle a_i^\dagger a_j^\dagger a_l a_k\rangle
\]

obeys (c^\dagger D c=\langle B^\dagger B\rangle\ge0), the particle-particle (two-particle) moment matrix (often called (P)). The hole analogue (B=\sum c_{ij}a_i^\dagger a_j^\dagger) gives (Q\succeq0). Both use degree-two operators and degree-four squares. But cubic (B=\sum c_{ijk}a_i a_j a_k) yields a 3-particle moment matrix (T), with six-operator entries: it requires a 3-RDM, not a 2-RDM cut. A symbolic CAR reducer should return a nonzero 3-body block for generic (c). A 2-RDM cut is valid only when higher-body coefficients cancel exactly or a lifted 3-RDM is retained.

## 4. The proposed numeric residual is not automatically a lower-bound certificate

If \hat A−bI = S+R with (S\succeq0) and ∥R∥ₒ ≤ η, then \hat A−(b−η)I \succeq0. The bound must be shifted by the operator-norm error; an unqualified “tiny residual” can reverse the inequality. A norm of coefficient residuals must be converted to a bound on the operator residual. Floating-point equality is also not a proof.

**Repair.** Use rational/algebraic coefficients with interval or exact directed-rounding arithmetic. Verify (R=R^\dagger) and bound ∥R∥ₒ rigorously. The operator Frobenius norm is itself a valid spectral upper bound; a coefficient-array Frobenius norm is valid only after accounting for the operator-basis Gram matrix. Set (b_{\rm safe}=b-η), and propagate solver tolerances into (E_{\rm lower}). Lean checks the finite identity but cannot make an uncertified norm estimate exact.

## 5. Rayleigh–Ritz gives an upper bound; residual gives neither ground-state identity nor a sharp gap

For a normalized trial state, (E_\psi=\langle H\rangle\ge E_0). A small residual ∥(H−Eψ)|ψ⟩∥ only proves concentration near eigenvalues; without a known spectral gap or an exclusion interval, the state may approximate an excited state or a mixture of nearby eigenstates. The inequality (E_{\rm lower}\le E_0\le E_{\rm upper}) remains sound when each side is independently certified; residual-based claims about the ground state do not.

**Repair.** Use residuals as diagnostics. Ground-state identification from an approximate eigenstate additionally needs justified spectral isolation; a positivity check only on its orthogonal complement does not remove off-diagonal couplings. The independently certified energy sandwich avoids this issue for energy accuracy. For a ground-state energy computation, the sandwich gap itself is sufficient; for observables, certify state fidelity or observable error using a gap-dependent theorem.

## 6. Fixed (N) is mathematically clean but chemically incomplete

The contraction from a 2-RDM to a 1-RDM uses fixed (N\ge2). Charge transfer within a molecule can remain in fixed total (N); grand-canonical systems and thermal ensembles may require multiple sectors. Spin, point-group, time-reversal, and spatial symmetries alter the affine constraints.

**Repair.** Keep fixed total (N) where appropriate, resolve relevant ((N,S,S_z)) and spatial-symmetry blocks, and for grand-canonical work solve (H-μ\hat N) on the direct sum or take the certified minimum over sectors.

## 7. A 2-RDM does not close dynamics or thermal physics

Two-body Hamiltonians make the *energy* linear in the 2-RDM, not the von Neumann entropy, finite-temperature free energy, or time derivative of the 2-RDM. The BBGKY equation couples the 2-RDM to the 3-RDM. Therefore the proposed representation is exact for ground-state optimization of two-body Hamiltonians, while reaction kinetics and finite-temperature quantities need an enlarged hierarchy or a controlled closure.

**Repair.** For dynamics, propagate a certified (k)-RDM hierarchy with a positivity-preserving closure and a posteriori residual bounds. For temperature, optimize over a consistent hierarchy including entropy bounds, or use thermofield/tensor-network primal witnesses and certify only energy/free-energy intervals that the dual actually proves.

## Falsifiable attack program

Use exact diagonalization as the oracle on a ladder of finite systems (Hubbard dimers, H chains, stretched geometries, frustrated plaquettes). For each Hamiltonian direction, solve SDP relaxations at increasing (N)-representability degree; ask the learned separator for an ((r,d)) certificate, then independently verify the CAR identity and norm margin with exact arithmetic. Measure:

\[
  \text{gap}=E_{\rm upper}-(E_{\rm lower}-\text{all certified tolerances}),
\quad (r,d,\text{bits},\text{support},\text{runtime}).
\]

The central hypothesis is falsified if, at fixed chemical accuracy, certificate degree or bit length grows rapidly with system size, or if adversarial directions (frustration, bond breaking, near-degeneracy) require the same rank/degree as the generic hierarchy. It is supported only if the sandwich gap closes with subexponential certificate resources while the proof survives independent exact verification.

The sharpest near-term test is a blind split: train the separator on small molecules, then freeze it and challenge it on stretched Hubbard ladders and random two-body directions. Success requires a valid lower bound, not merely an accurate energy prediction.

## Appendix: an exact extensive DQG gap diagnostic

The proposed disjoint-triple family is mathematically sound as an elementary stress test, subject to explicitly checking the full chosen (D,Q,G) index conventions. Take (2r) disjoint triples of spinless modes, (N=3r), and (H=\sum_T\sum_{i<j\in T}n_i n_j). Set the diagonal pseudo-marginal (p_i=1/2), (q_{ij}=1/8) within a triple and (1/4) across triples, with off-diagonal 2-RDM entries zero. Its fixed-(N) contraction is

\[
\sum_{j\ne i}q_{ij}=2(1/8)+(6r-3)(1/4)=(3r-1)/2=(N-1)p_i.
\]

The (D/P) and (Q) diagonal entries are nonnegative. The number-operator block of the particle-hole G matrix can be written

\[
P=\tfrac14J+\operatorname{blockdiag}\!\left[\tfrac38I_3-\tfrac18J_3\right]\succeq0,
\]

since each triple block has eigenvalues (0,3/8,3/8), and the global (J) contributes only in the all-ones direction. The candidate energy is (2r\cdot3\cdot(1/8)=3r/4). The exact CAR identity

\[
C_T=\big\{a_{T,1}a_{T,2}a_{T,3},(a_{T,1}a_{T,2}a_{T,3})^\dagger\big\}=1-S_T+\sum_{i<j\in T}n_i n_j
\]

gives \sum_T C_T=H-r on the (N=3r) sector, hence every physical state has (H\ge r). Separately, \sum_T(S_T-3/2)^2=2H-3r/2\ge0) gives the DQG-compatible lower bound (H\ge3r/4), exactly matching the pseudo-marginal objective. Thus an extensive (r/4) DQG gap survives while a linear-size cubic anticommutator certificate closes it.

The root synthesis and executable diagnostic now expand (G_{ij,kl}=\langle a_i^\dagger a_j a_l^\dagger a_k\rangle) in an explicitly stated Gram convention, checking every entry for six modes and providing the block formulas for arbitrary r. The full G block audit passes with exact rational arithmetic; see marginal_hunt_synthesis.md and ../results/marginal_hunt_witness.json. This is an elementary diagnostic family, not evidence about generic chemistry.
