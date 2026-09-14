# Intervention-stable Pauli representations

## Verified primary sources

- Hazime Mori, “Transport, Collective Motion, and Brownian Motion,” *Progress of Theoretical Physics* 33 (1965), 423–455, DOI `10.1143/PTP.33.423`: https://doi.org/10.1143/PTP.33.423
- Robert Zwanzig, “Memory Effects in Irreversible Thermodynamics,” *Physical Review* 124 (1961), 983–992, DOI `10.1103/PhysRev.124.983`: https://doi.org/10.1103/PhysRev.124.983
- E. Bairey, I. Arad, and N. H. Lindner, “Learning a Local Hamiltonian from Local Measurements,” *Physical Review Letters* 122, 020504 (2019), DOI `10.1103/PhysRevLett.122.020504`: https://doi.org/10.1103/PhysRevLett.122.020504
- Riccardo Massidda, Atticus Geiger, Thomas Icard, and Davide Bacciu, “Causal Abstraction with Soft Interventions,” *Proceedings of the Second Conference on Causal Learning and Reasoning*, PMLR 213:68–87 (2023): https://proceedings.mlr.press/v213/massidda23a.html

## Bounded result and experiment

For Pauli observables (O_i), an exactly closed reduced model exists when their span is invariant under every controlled adjoint generator \(\mathcal L_u^\dagger\). Otherwise write

\[
\mathcal L_u^\dagger O_i=\sum_j A_{ij}(u)O_j+R_i(u).
\]

The implementation computes exact rational Pauli coefficients of (i[H,O]), the projected generator matrix, and the omitted-support (l_1) residual. If ‖A‖ is bounded by (L), a Euclidean comparison gives

\[
\|z(t)-\hat z(t)\|_2\le e^{Lt}\delta_0+\frac{e^{Lt}-1}{L}\sup_{s\le t}\|r(s)\|_2.
\]

In a normalized Pauli-orthogonal basis, the row generator is skew-symmetric, (A^T=-A). Its propagator has Euclidean norm one, so the error bound improves to

\[
\|z(t)-\hat z(t)\|_2\leq\delta_0+\int_0^t\|r(s)\|_2\,ds
\leq\delta_0+t\epsilon,
\qquad \epsilon=\sup_s\sqrt{\sum_i b_i(s)^2},
\]

Here $b_i=\sum_{P\notin V}|r_{iP}|$ is the sum of the absolute omitted Pauli coefficients in row $i$, which bounds $|\operatorname{Tr}(\rho R_i)|$ uniformly over all density matrices. Individual coefficients alone would not bound the whole row. If the retained span is invariant under **each** control generator, it remains invariant under every piecewise bounded time-dependent control, since each interval evolves within the same span and the propagator products preserve it.

A concrete counterexample is (H_A=XII) with retained observables ({ZII,YII}): (i[H_A,ZII]=2YII) and (i[H_A,YII]=-2ZII), so this two-dimensional space closes. Adding the unseen coupling (H_B=ZXI) gives

\[
i[H_B,YII]=2XXI,
\]

which is omitted from the retained space and falsifies intervention-stable closure immediately.

A minimal validator uses exact dense matrices for (n\le6). A valid representation must retain its predeclared residual bound under a withheld coupling or expose a missing observable/memory term. These closure and error statements are elementary consequences of finite-dimensional linear dynamics; no new theorem is claimed here.

## Limits

Mori–Zwanzig theory makes memory/noise from discarded variables explicit; it does not imply a compact Markovian closure. Hamiltonian learning guarantees apply to declared locality/model classes and do not establish universal intervention-stable abstraction. A latent predictor that fails under a new coupling is not a causal abstraction certificate.
