# Ground-energy certificates from exact excitation counting

Two connected H6 Hamiltonians now admit tighter ground-energy intervals with fewer response directions. The proof certifies a coarse second-eigenvalue threshold in the exactly spin-symmetric, Sz=0 sector and turns the physical witness's small variance into a sharp ground lower bound. Original-H intervals remain limited by the separately checked5.4e-11-Hartree spin-symmetry perturbation allowance.

| Perturbation | Previous Schur directions | Temple prefix directions | Hs interval width (Ha) | Original-H interval width (Ha) |
|---|---:|---:|---:|---:|
| t=1/100, coefficient norm |31|27|1.7372603689662484e-14|5.401737260369442e-11|
| t=1/50, square norm |31|28|1.2468951846429005e-13|5.412468951848323e-11|

A subsequent fresh t=1/50 construction also passes, using29 newly proposed response directions and the same final interval. It imports the direct reference P/factors and initial upper seed, then independently refines the upper and generates its response from zero directions. It does not import the31-direction energy proof or the28-direction compressed prefix.

These prefix results inherit their response vectors and upper witnesses from accepted connected transfers. Prefix compression took24.56 and26.30 seconds, including exact acceptance. It does not establish cheaper initial discovery. Certificates occupy2,228,784 and2,247,822 bytes; the explicit reference factors dominate their size. Both replays still evaluate400 distinct configurations. Fewer response directions do not mean a compressed many-body state space.

## The eigenvalue-count argument

All eigenvalues in this argument belong to Hs restricted to Sz=0. The existing complete reference factors and exact perturbation norm establish

\[
QH_sQ\succeq\gamma Q,\qquad \gamma>\tau.
\]

The residual-response machinery constructs a full-P matrix S_K(tau) satisfying

\[
S_K(\tau)\preceq S(\tau)
=P(H_s-\tau)P-PH_sQ[Q(H_s-\tau)Q]^{-1}QH_sP.
\]

This is the same global Loewner lower bound used by the [response certificate](marginal_sparse_response.md); it is not a condition checked only on the selected response vectors. The response directions improve the bound, while S_K still acts on every retained P coordinate.

Strict Q positivity permits Schur congruence. Loewner monotonicity then implies

\[
n_-(H_s-\tau)=n_-(S(\tau))\le n_-(S_K(\tau)).
\]

If the final exact count is at most one, an independently evaluated physical witness with Rayleigh value mu<tau forces exactly one eigenvalue below tau. Consequently lambda2>=tau. Zero eigenvalues at tau are allowed; they must be counted as zero, not as negative or positive.

The checker computes inertia using rational congruences. It eliminates any nonzero diagonal as a1x1 block. If all remaining diagonals vanish but an off-diagonal entry b is nonzero, the block [[0,b],[b,0]] contributes one positive and one negative eigenvalue. A remaining zero matrix contributes its full dimension to the nullity. Floating eigenvalues select proposals only.

## Temple lower bound and Hamiltonian transfer

For a normalized physical witness psi, replay recomputes both

\[
\mu=\langle\psi,H_s\psi\rangle,\qquad
v=\|(H_s-\mu)\psi\|^2.
\]

Since the spectrum has exactly one eigenvalue E0 below tau, (Hs-E0)(Hs-tau) is positive semidefinite. Taking its expectation yields

\[
v\ge(\mu-E_0)(\tau-\mu),\qquad
\boxed{E_0(H_s)\ge\mu-\frac{v}{\tau-\mu}}.
\]

This is the established Temple inequality, not a new eigenvalue theorem; see [Cape, Tang and Priebe](https://arxiv.org/abs/1603.06100). The implemented contribution combines exact residual-response inertia, checked molecular reference factors, physical variance, and independent rational replay. A loose supplied upper endpoint is never substituted for the actual Rayleigh value in the numerator.

Exact SU(2) symmetry and even particle number transfer the lowest Sz=0 energy to the global fixed-N minimum. They do **not** transfer the Sz=0 second eigenvalue or nondegeneracy to the full fixed-N sector; spin multiplets can create repeated global ground energies.

Finally, replay recomputes epsilon>=||H-Hs|| from exact CAR coefficient differences. It exports

\[
\mu-\frac{v}{\tau-\mu}-\epsilon
\le E_0(H)\le\langle\psi,H\psi\rangle.
\]

The upper is evaluated directly on H. The reference-to-Hs gap norm and Hs-to-H energy allowance remain separate.

## Exact finite results

At t=1/100, tau=-1266420690003/200000000000. The27-direction Schur lower matrix has exact inertia31 positive,1 negative,0 zero. The preceding26-direction prefix has two numerical negatives and is not accepted as an excitation certificate. The witness variance is1.7372603679660957e-17 Ha².

At t=1/50, tau=-3166119057599/500000000000. The28-direction Schur lower matrix again has exact inertia(31,1,0). Its witness variance is1.2468951845353828e-16 Ha². The exact CAR square norm is required to justify the retained reference threshold; deleting that square certificate causes replay to reject the more optimistic threshold under the default coefficient-sum norm.

Both thresholds lie approximately0.001 Ha above their exact witness Rayleigh values. These are certified lower bounds on the second Sz=0 eigenvalue, not computed second eigenvalues. Both energy certificates have independently passed standard-library-only replay:

```sh
python -S -m experiments.marginal_spin_temple \
  --verify results/marginal_h6/connected_spin_1_100_temple/certificate.json
python -S -m experiments.marginal_spin_temple \
  --verify results/marginal_h6/connected_spin_squares_1_50_temple/certificate.json
```

Tests cover singular inertia matrices,2x2 pivots, known inertia under invertible congruence, threshold degeneracy, zero variance, witness normalization and support outside P, mu approaching tau, false excitation counts, invalid symmetry/sectors, missing Q coverage, separate original-H transfer, and the combined reference-norm proof.

## Remaining construction problem

A fresh bounded proposer now targets the second Schur eigenvector at a fixed tau, uses sparse Q actions for the inverse-response proposal, and stops only after exact inertia and full interval replay succeed. The second-vector selection is a heuristic. Two orthogonal vectors can both have negative quadratic forms even when the true matrix has only one negative eigenvalue, so this selection rule alone is not a convergence proof. Response-cap failure, insufficient variance accuracy, or loss of the strict complementary gap must produce no accepted certificate.

The t=1/50 fresh run starts with zero response directions and finishes at29. Its recorded transfer phase takes125.95 seconds, versus404.28 seconds for the earlier full-energy Schur response on this perturbation. These are individual observed runs with differing stopping conditions, not a controlled performance benchmark. The fresh path avoids repeated full-energy bisection and intermediate factor replays, while final acceptance recomputes every premise. Final replay still evaluates368 reference-H0,400 current-Hs, and400 original-H source states drawn from400 configurations, and the reference factors remain200 and168 dimensional. No factor of the joined368-state Q block is constructed.

The fresh final Schur lower matrix has exact inertia(31,1,0), even though its numerically smallest eigenvalue is about-25.8563. That large negative direction is allowed: the theorem needs an eigenvalue count, not an accurate ground-energy estimate from the Schur matrix itself. The next numerical eigenvalue is about0.00072437; neither floating value is used for acceptance. The precise lower endpoint comes from the exact physical variance and the certified separator. This separates coarse global spectral exclusion from fine ground-energy accuracy. The fresh certificate occupies2,261,572 bytes, so the reduced directional work has not removed the dominant reference-factor storage.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_spin_transfer \
  --source results/marginal_h6/direct_spin/certificate.json \
  --output results/marginal_h6/connected_spin_squares_1_50_temple_fresh \
  --strength=1/50 --norm-squares --temple-offset=1/1000

python -S -m experiments.marginal_spin_temple \
  --verify results/marginal_h6/connected_spin_squares_1_50_temple_fresh/certificate.json
```

The output directory must be new when rerunning construction. A focused integration test makes the old full-energy response constructor and new-factor proposer raise if called; the fresh Temple transfer still passes. A deliberately inaccurate upper fails with a saved rejection and no final certificate.

The unresolved large-system problem remains obtaining useful reference gaps, physical witnesses, and exact response data without explicit growth in configuration support. This finite excitation-count route reduces the demanded spectral accuracy of the lower proof; it does not solve universal physical-marginal representability.
