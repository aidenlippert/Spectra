# Exact singlet energy integration for the eight-site Hubbard chain

Later work tightens the interval to [-4.23585,-4.23580593] using exact symmetry-orbit moments, with a slightly sharper range-three variant. See [the current moment report](marginal_symmetry_moments.md). The following records the preceding determinant-response milestone.

The fixed eight-site half-filled open Hubbard model at U=4 and t=1 now has an accepted exact ground-energy interval:

\[
-4.3\le E_0\le
-\frac{525480580781980877341375}{125000000011284578672289}
<-4.20384.
\]

The exact width is approximately 0.0961553541 in hopping-energy units. This closes the finite-model gap between the compact charge-complement certificate and an actual energy interval. It is not a general quantum-chemistry solver or a chemical-accuracy result.

## Representation and proof gates

Fourteen noncrossing valence-bond vectors replace the 70 determinant coordinates of the balanced valence sector. The independent embedding verifier reconstructs their coefficients, computes the CAR spin-raising action, checks S+V=0, computes rank(S+)=56, and verifies a positive 14-by-14 Gram matrix. Thus the vectors span the entire valence singlet subspace. Six mutations, including changed coefficients and a duplicated basis vector, are actually submitted to the verifier and refused.

The basis is nonorthogonal. Every Schur calculation retains its exact Gram G=VᵀV; it never substitutes an identity matrix or treats these as 14 determinant states.

The ground-spin step uses [Lieb's repulsive half-filled bipartite Hubbard theorem](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.1201). The replay requires the original Hamiltonian, positive U=4, real spin-independent nearest-neighbor hopping t=1, a connected open chain, half filling, and balanced sublattices of four sites each. These hypotheses imply a singlet ground state. The established theorem is an explicit mathematical dependency, not reproved by this artifact; SU(2) commutation alone is not used to infer ground spin.

Within S=0, the orthogonal complement of the complete valence singlet space is entirely ionic (D>=1). The previously accepted charge-product DP certificate therefore supplies A=QHQ >= gamma=-3.81 on this complement. The energy replay reruns that DP certificate and verifies the same original Hamiltonian. It also checks exact spin commutators.

Write W=QHV. The response consists of 14 columns B=p_5(A)W, where p_5 has six rational coefficients. Polynomial application preserves the exact singlet space and applies Q at every step. Coefficients are discovered numerically by fitting an inverse approximation in Chebyshev coordinates. No numerical fit or spectrum is accepted as proof.

For a proposed lower endpoint tau<gamma, replay forms

\[
D=B^\top[A^2-(\tau+\gamma)A+\tau\gamma I]B,
\qquad L=B^\top(A-\gamma I)W,
\]

and verifies D>0, an exact solve, and exact positive LDL pivots for

\[
S=-\tau G-\frac{W^\top W-L^\top D^{-1}L}{\gamma-\tau}.
\]

This follows from the resolvent upper bound

\[
W^\top(A-\tau I)^{-1}W
\preceq\frac{W^\top W-L^\top D^{-1}L}{\gamma-\tau},
\]

which is valid for A>=gamma. Positive S proves the lower endpoint. An independently supplied integer vector is evaluated under the original CAR Hamiltonian for the variational upper endpoint.

## What was resolved, and what was not

The previous rank>=38 obstruction concerned the full 70-coordinate valence reference and a fixed target. It does not apply after restricting to the proved ground-spin sector and changing the retained geometry. The new representation uses 14 retained vectors and 14 response vectors without raising the existing oracle's 4096-action cap.

The upper witness also proves that tau=-4.2 cannot be a ground-energy lower endpoint: an explicit physical state has energy below -4.20384. Earlier rank diagnostics at -4.2 remain valid diagnostics, but that target is now ruled out independently of response quality.

The accepted replay still takes 3870 distinct determinant source actions and references 4234 determinants, out of 4900 in the balanced-spin sector. Each response column has 3476–3716 nonzero coefficients. The short polynomial recipe does not yet produce a compact contraction cost. This is the immediate remaining representation problem: evaluate the required moments without expanding nearly the whole determinant sector.

The separate ten- and twelve-site charge-DP results establish weaker complement bounds only. This singlet energy implementation is deliberately fixed to the original eight-site model. Transfer to nonuniform molecular Hamiltonians, equal-accuracy scaling, thermal/response quantities, and general physical-marginal representability remain open.

## Artifacts and replay

- `results/marginal_graded_hubbard8/singlet_energy_integer/certificate.json`: accepted embedded gap/basis/polynomial/upper recipe.
- `results/marginal_graded_hubbard8/singlet_energy_integer/receipt.json`: exact construction acceptance.
- `results/marginal_graded_hubbard8/singlet_energy_integer/independent_replay.json`: separate standard-library replay output.
- `results/marginal_graded_hubbard8/discovery/singlet_energy_replay.py`: read-only exact research verifier.
- `results/marginal_graded_hubbard8/discovery/discover_singlet_energy.py`: numerical proposal driver; exact replay is mandatory.
- `results/marginal_graded_hubbard8/singlet_valence_embedding/root_independent_replay.json`: basis completeness and tamper rejection.
- `tests/test_marginal_singlet_h8_research.py`: four integration/refusal tests.

Run from the repository root:

```sh
python -S results/marginal_graded_hubbard8/discovery/singlet_energy_replay.py results/marginal_graded_hubbard8/singlet_energy_integer/certificate.json
python -S -m unittest tests.test_marginal_singlet_h8_research -v
```

The initial embedding exporter and its first agent-written replay did not establish all claimed checks. Their receipt is marked superseded. The current independent replay directly verifies spin annihilation, completeness, the actual H_PP block, and exact LDL reconstruction. An intermediate numerical response probe also omitted Q in one action used to form the squared moment; that omission was corrected before these energy certificates were constructed. No result from either incomplete verifier is used as acceptance evidence.

## Next contraction target, with an exact algebra check

The response need not carry explicit Q insertions through a tensor network. Define full moments M_n=VᵀHⁿV and normalized matrices A_n=G⁻¹M_n. Set C_0=I and

\[
C_n=-\sum_{j=1}^{n} A_j C_{n-j}.
\]

Then the projected moments are

\[
K_k=V^\top H Q(QHQ)^k QHV=-G C_{k+2}.
\]

For example, K_0=M_2-M_1G⁻¹M_1. This follows by expanding R(z)=Vᵀ(z-H)⁻¹V and **T(z)=zG⁻¹R(z)**, then applying the nonorthogonal Schur identity

\[
K(z)=zG-M_1-GR(z)^{-1}G.
\]

For a degree-five response, all required matrices reduce to K_0 through K_12, hence M_0 through M_14. An exact four-site CAR diagnostic independently computed both routes and found equality for every entry through K_12 (`projected_moment_recurrence/receipt.json`). That check validates the recurrence, not a fast moment-contraction implementation.

The next implementation target is therefore precise: contract the full moments M_0 through M_14 between prescribed noncrossing valence-bond boundaries, without expanding determinant vectors. A local-term expansion is polynomial in site count for fixed moment order **per prescribed pair of boundary vectors**, but its exponent grows with order. The complete noncrossing singlet basis has Catalan growth, and long nested bonds can also enlarge tensor bond dimensions. Neither this recurrence nor an MPO description proves efficient general scaling.


## Final validation

All 386 marginal tests passed in 279.336 seconds. Separate `python -S` replay accepted the final six-integer polynomial recipe, with coefficients `[241387155216,-60004966197,10665249314,-1043116424,49452130,-886207]` in ascending power order. Its recorded construction replay took 21.839 seconds in this CPU run. Scaling the response polynomial has no effect on its span; rounding the original rational coefficients was followed by complete exact acceptance. Both earlier rational and final integer certificates are preserved.
