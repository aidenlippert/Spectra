# Primary sources for interacting-certificate scaling (2026-09-15)

## Cluster and coarse-grained precedents

* Lin & Lindsey, “Variational embedding for quantum many-body problems,” [arXiv:1910.00560](https://arxiv.org/abs/1910.00560). The lower bound is the optimized SDP value (or a checked feasible dual certificate), not an arbitrary feasible relaxed marginal. The method uses local/global linear and PSD constraints for spin and fermionic algebras.
* Khoo & Lindsey, “Scalable semidefinite programming approach to variational embedding for quantum many-body problems,” [arXiv:2106.02682](https://arxiv.org/abs/2106.02682). Clusters are glued by explicit global consistency updates; local PSD blocks and overlap maps \(\mathcal R_{C\to O}(X_C)=X_O\) are the reusable pattern.
* Kull, Schuch, Dive & Navascués, “Lower Bounding Ground-State Energies of Local Hamiltonians Through the Renormalization Group,” [arXiv:2212.03014](https://arxiv.org/abs/2212.03014). Coarse-graining maps remove many consistency constraints while preserving a convex relaxation; linear optimization over that relaxation yields the bound. Tightness depends on the chosen scheme.
* Li & Lu, “Quantum variational embedding for ground-state energy problems: sum of squares and cluster selection,” [arXiv:2305.18571](https://arxiv.org/abs/2305.18571). Gives an SOS SDP hierarchy and cluster-selection strategies; selection improves tightness heuristically, while exact PSD/consistency constraints carry validity.

For \(M=\sum_aB_a^\dagger G_aB_a\), \(G_a\succeq0\), each term is PSD and its nullspace is \(\ker(\sqrt{G_a}B_a)\), not generally \(\ker(B_a)\) when \(G_a\) is singular. Under \(B'_a=T_aB_a\), use \(G'_a=T_a^{-\dagger}G_aT_a^{-1}\); congruence preserves positivity, while a non-invertible map may change the cone.

## Spin trace

“Heisenberg models and Schur--Weyl duality,” [arXiv:2201.10209](https://arxiv.org/abs/2201.10209), supplies the SU(2) decomposition. For even (N),
\[
\mathcal H_N=\bigoplus_S V_S\otimes\mathcal M_S,
\qquad A=\bigoplus_S I_{V_S}\otimes A_S
\]
for spin-invariant (A). Since (M_s=0) contains one weight per integer (S), and (M_s=1) one per (S\ge1),
\[
\operatorname{Tr}_{N,S=0}A=\operatorname{Tr}_{N,M_s=0}A-\operatorname{Tr}_{N,M_s=1}A.
\]
For general (A), use \(\bar A=\int_{SU(2)}U A U^\dagger dU\); (P_0) commutes with (U), so \(\operatorname{Tr}(P_0A)=\operatorname{Tr}(P_0\bar A)\). Twirling changes nonsinglet components and must be included in the certified operator.

## Existing acceptance equations

`experiments/marginal_symbolic.py:verified_residual` constructs
\[
R=H-bI-\sum_aQ_a^\dagger Q_a-(N-N_0)X,
\qquad \eta=\sum_w|R_w|,
\]
and returns \(b-\eta\), using normalized CAR moments bounded by one. `spin_replay.check` twirls (H) and charges \(\delta=\|H-\operatorname{twirl}(H)\|_{1,\mathrm{coeff}}\).

`spin_screen.check_sector` uses \(Z=\alpha_shift(m,n)-M_s\) and the code identity
\[
H_{\rm twirl}-(ZY)-aS^2-[S_+W+(S_+W)^\dagger]
=bI+\sum_aQ_a^\dagger Q_a+(N-N_0)X+R.
\]
The ladder term is `ladder_ideal(W)=S_+W+(S_+W)^dagger`, not an anticommutator \(\{L^\dagger,L\}\). Casimir/ladder terms are sector ideals. `spin_screen.check` combines singlet and (M_s=1) receipts by the minimum, then charges the spin defect once.

`commutator_dual_witness.check` requires \(y_{()}=1\), \(y(w)=y(w^\dagger)\), \(|y_w|\le1\), all \(y((N-N_0)p)=0\) for `multiplier_basis(max_body=2)`, and exact PSD Grams \(G_{ij}=y(Q_i^\dagger Q_j)\succeq0\). `propose` performs rational affine elimination and trace-seed mixing; every candidate must pass `check`. This is a diagnostic, not a numerical certification.

A bounded direct singlet-dual repair should parameterize \(y_w\), impose normalization, Hermiticity, exact singlet number/Casimir/ladder ideal equalities, and coefficient boxes before trace mixing; build twirled local Grams and require exact rational PSD; enforce forced-null rows/columns before mixing with the singlet trace seed. Report objective and refusals only until a complete exact replay accepts every constraint.
