# Stretched H4 square: accepted molecular certificate

The square H4 fixture has side length 2.0 Angstrom, uses STO-3G RHF canonical orbitals, and contains 8 spin orbitals and 4 electrons. PySCF FCI gives electronic energy -3.330388605150463 Ha. An independent 70-state CAR diagonalization of the rational exported Hamiltonian gives -3.330388605150749 Ha. This agreement checks the encoding numerically; it does not certify the underlying integral algorithm or basis approximation.

The exact rational certificate now brackets its electronic ground energy:

```text
-3.3303887263976906 <= E0 <= -3.330388605150749
width = 1.212469417001217e-7 Ha
```

The lower endpoint comes from rational CAR square expansion and residual accounting. The upper endpoint is independently recomputed from an integer-amplitude wavefunction. The displayed upper decimal is rounded; the receipt stores the exact rational endpoint.

The unsplit spin-charge formulation was interrupted after more than ten minutes in Clarabel without an accepted proposal. It is recorded as an interrupted computation, not infeasibility. An exact scan of Hamiltonian monomial support found seven nonidentity binary symmetry masks (three independent bits, including alpha/beta parity). Splitting the words by those characters and filtering the coefficient equations reduced the largest Gram block from 92 to 46. The accepted formulation has 50 blocks, 435 coefficient equations, and 111 multiplier directions. Assembly took 0.57 seconds and the numerical solve 6.94 seconds. Exact replay is the acceptance criterion, regardless of solver status.

The rational coefficient denominator is 10^12. The saved fixture's conservative coefficient perturbation enclosure is 4.631027315316937e-10 Ha. It includes all raw coefficients potentially dropped by the original 1e-13 threshold, plus exact retained-coefficient rounding error. This error concerns export of computed floating integrals; integral algorithm, geometry, finite-basis, and physical-model errors are not certified. The generator has been corrected to retain all raw coefficients in future error calculations.

Reproduce the proposal from the saved fixture:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_molecule_stress --solve-certificate
```

Replay the accepted interval without numerical dependencies:

```sh
python3 -S -m experiments.marginal_transfer_verify results/marginal_molecule_stress/accepted_degree3_certificate.json
```

This is a second tested molecular geometry, not evidence of general molecular scaling. The measured run initially scanned 256 binary masks. The generator now obtains the same three independent symmetries by GF(2) elimination of the monomial-support equations, avoiding exhaustive mask enumeration. A test compares its full generated group against the original exhaustive eight-mode scan.
