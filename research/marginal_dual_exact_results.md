# An exact nonrepresentable cubic moment witness

The ten-mode, five-particle matched Hamiltonian at t=1/5 now has a rational moment witness proving that its degree-six moment relaxation is strictly looser than the physical problem.

The witness has energy

\[
L(H)=\frac{29314712562242917}{9000000000000000}
\approx3.257190284693658.
\]

The independent exact sector calculation established

\[
E_0>\frac{1803752506297}{549755813888}
\approx3.281006695573524.
\]

Thus the relaxation gap is strictly greater than 0.0238164108798. This conclusion does not depend on a numerical solver claiming optimality. The witness need not be optimal: exact feasibility and its energy below the physical ground energy suffice.

## Which relaxation is certified

The functional L is defined on all CAR polynomials through ladder degree six. It is normalized, Hermitian, invariant under the independent matched-pair occupation phases, and satisfies every fixed-number relation (Nhat-5)X with real Hermitian X through ladder degree four. Its moment matrix is positive semidefinite on every operator polynomial B through degree three.

To implement the last condition, generate every normal-ordered monomial of degree zero through three, including overlapping creation/annihilation indices, and group by its full conserved charge vector. The resulting 221 blocks cover the entire operator space through degree three. Positivity of each block implies positivity of their block-diagonal direct sum, including complex linear combinations. This is stronger than trusting a certificate-supplied list of favorable blocks; the verifier regenerates the complete list.

Any exact SOS identity H-b=sum B†B+(Nhat-5)X at these degree limits obeys L(H)>=b. Consequently no exact certificate in this relaxation can reach the physical ground energy. The witness alone supplies the needed upper bound on this exact SOS objective; no strong-duality assumption or claim about the relaxation's optimum is needed. Existing residual-corrected certificates remain valid physical energy bounds and are not used here to assert a matching relaxation optimum.

The claim applies to this explicit instance and these degree limits. It does not prove a universal particle-number transition, a universal quartic insufficiency result, or optimality of the displayed witness. In particular, a formulation with higher-degree number multipliers must be treated separately.

## Exact reconstruction

Clarabel first proposes a numerical equality-dual vector from the coefficient-space SOS problem. The 541 retained coordinates correspond to one representative of each invariant Hermitian monomial pair. An off-diagonal representative coordinate is twice the functional value of either member; diagonal representatives use the coordinate directly. The verifier implements this convention explicitly.

The proposer builds the rational affine system consisting of normalization and the number relations. It computes exact reduced row-echelon form, rounds only free coordinates to denominator 10^10, and solves pivot coordinates exactly. This preserves the forced moment-kernel identities.

It then mixes the reconstructed vector with the maximally mixed state on the full fixed-N sector by weight 1/100000. For canonical diagonal rank-k words, that physical state's moments are

\[
(-1)^{k(k-1)/2}\frac{\binom Nk}{\binom Mk};
\]

off-diagonal moments vanish. The mixture supplies a small positive margin on the nonmandatory supports while preserving the affine equalities. It raises the energy slightly but leaves a large separation from the physical minimum.

The final acceptance check uses only rational arithmetic. Every moment block undergoes exact symmetric LDL elimination. Positive pivots contribute their Schur complement; a zero pivot is accepted only when the remaining column is exactly zero; a negative pivot rejects the witness. No floating eigenvalue tolerance appears in acceptance.

## Verification receipt

- 541 rational coordinates.
- 1,091 generated Hermitian number-multiplier basis relations checked, including relations that vanish by symmetry.
- 221 generated moment blocks proved positive semidefinite exactly.
- Matched Hamiltonian and particle number checked against the declared problem.
- Standard-library-only replay succeeds under `python3 -S`.

Tests compare physical reference moments with independent signed fermionic traces, exercise singular and indefinite PSD cases, reject normalization/Hamiltonian tampering, and compare the saved witness energy against the independent spectral lower bound. The numerical helper also tests its equality-dual convention against an exactly soluble number-operator problem.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_dual_exact --build
python3 -S -m experiments.marginal_dual_exact --verify results/marginal_dual/witness.json
```

Implementation: `experiments/marginal_dual_exact.py` and `experiments/marginal_adaptive.py`. The exact artifact is `results/marginal_dual/witness.json`; `receipt.json` records the check, and `numeric.json` records the initial numerical diagnostics.
