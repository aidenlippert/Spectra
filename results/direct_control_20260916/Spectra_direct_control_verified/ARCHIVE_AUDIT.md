# Audit of the uploaded certified-control bundle

## Finding

The uploaded bundle is a valid-looking finite-dimensional control certificate for the declared H8 instance, subject to its stated assumptions and the parent agent's fresh replay. Its Hamiltonian is physically equivalent to the earlier canonical H8 Hamiltonian by spatial-orbital phase flips. It is not evidence of a new physical molecule or a general many-body solution.

The original inputs are canonical RHF H8 data: eight hydrogen atoms at 1.4 Angstrom spacing, STO-3G, 16 spin orbitals, eight electrons, and the balanced \((N_\alpha,N_\beta)=(4,4)\) sector. The uploaded fixture's `orbital_basis` is `RHF canonical MO`; its nested `localized_basis` metadata says localization was not run.

## Exact Hamiltonian equivalence

The uploaded and earlier `h8_cold` coefficient files use different serialized rational values, but a direct 256-mask test over spatial orbital signs found the exact equivalence

\[
a_{p\sigma}\mapsto s_p a_{p\sigma},
\qquad
(s_0,\ldots,s_7)=(1,-1,1,1,1,-1,1,1).
\]

For every common Hamiltonian word, this mask maps the earlier coefficient to the uploaded coefficient exactly. The maximum and total absolute coefficient residuals are both zero at the declared rational precision. The earlier coefficient mismatch therefore reflects an orbital phase convention, not a demonstrated difference in the underlying physical Hamiltonian.

This does not by itself establish equivalence of the complete driven certificate. The phase transformation leaves \(D=n_6+n_7-n_{10}-n_{11}\) unchanged, but it changes the sign of the hopping terms in \(W\) involving spatial orbital 5. The MPS coefficients and the control waveform must be transformed consistently before claiming identical trajectories or identical control certificates. The correct claim is an exactly phase-equivalent Hamiltonian instance with a separately supplied state/control convention.

## Defect certificate

For a piecewise differentiable candidate path \(q(t)\), with \(r(t)=i\dot q(t)-K(t)q(t)\), Hermitian unitary evolution gives the error bound from the initial mismatch, the integral of \(\|r(t)\|\), and the jumps between segments. No exponential factor involving \(\|H\|\) is needed. The implementation reconstructs the physical Hamiltonian and control actions from rational inputs, forms exact small Gram matrices, and evaluates continuous polynomial residual integrals with outward rounding. It does not accept sampled residuals or an assumed Galerkin identity. The observable bound correctly keeps the candidate's final norm explicitly, so a non-unit-norm \(q\) is handled rather than silently normalized.

The accepting code checks the fixture and MPS hashes, Hermiticity, fixed spin populations, MPS charge flow, control norm bounds, reduced-matrix symmetry, dimensions, positive durations, declared control points, and malformed trajectories. The eight-coordinate case is rejected because its robust interval does not prove the target. That is an insufficient-proof result, not a physical impossibility result.

## Robustness allowance

For the declared model, the combined observable allowance is

\[
16Ta+4r_0+4\Gamma,
\]

where \(a\) bounds each control-coefficient perturbation pointwise, \(r_0\) is standard trace distance to the normalized supplied MPS, and \(\Gamma=\int_0^T\sum_j\gamma_j(t)\,dt\) for the specific Lindblad terms \(\gamma_j(Z_j\rho Z_j-\rho)\). Since \(\|D\|\le2\), waveform uncertainty contributes \(16Ta\); initialization contributes \(4r_0\); and the phase-flip channel contributes \(4\Gamma\). Thus the short-case allowance \(19/500\) and long-case allowance \(453/12500\) are arithmetically consistent with their declared budgets.

These are conditional mathematical guarantees. The bundle does not measure preparation fidelity, calibrate a laboratory pulse, or cover arbitrary environmental noise. A mixed initial state is covered only when it remains in the same balanced-spin sector and lies within the declared trace-distance ball around the supplied MPS.

## Sector and cost scope

The checker reconstructs the exact action of the supplied Hamiltonian on the balanced 4-alpha/4-beta sector. It does not construct the full eight-electron Fock-space matrix or certify other spin-population sectors. “All original terms remain present” means all supplied terms are retained in this accepted sector action.

The parent agent's fresh original replay completed in 47.8699 seconds with 173,621,248 bytes peak memory, with the expected accepted cases and deliberate refusal. Discovery and acceptance enumerate the 4,900-state balanced sector and construct the full sector action; the long basis contains 313,600 integer coefficients. The small online Gram matrices therefore do not establish cheap discovery or scalable many-body compression.

## Independent next discovery test

For the short \(T=2\) pulse, freeze the supplied MPS, controls, four segment boundaries, and observable. Perform a fresh residual-Krylov discovery beginning from the MPS and applying \(H,D,W\) at the segment centers. Orthogonalize and rank added directions by predicted reduction of the short-pulse defect under a fixed action budget. Do not import the supplied 24-dimensional embedding, its trajectory coefficients, or its full-trajectory snapshots. Let the rational checker decide acceptance, and use full-sector propagation only as a post hoc numerical comparison. Compare any accepted construction to the 24-coordinate result at matched discovery cost.
