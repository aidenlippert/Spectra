# Transfer after breaking flavor permutation symmetry

A compact proof now supports unequal pair hoppings by retaining its 21 averaged square orbits and adding eight unaveraged linear squares. Both tested Hamiltonians have exact rational lower and upper bounds. The old source proof is not expanded into thousands of individual images.

This is **analytic perturbation transfer**, not a new discovery algorithm for asymmetric Hamiltonians. Its intervals are substantially wider than the matched source interval. That loss of tightness is the unresolved problem exposed by the test.

## Hamiltonians and exact identity

Starting from the ten-mode, five-particle matched Hamiltonian, change the pair hoppings by

\[
\Delta H=-\sum_{i=0}^{4}e_i T_i,\qquad
 e_i=\epsilon(i-2),\qquad
 T_i=a_i^\dagger a_{i+5}+a_{i+5}^\dagger a_i.
\]

The five hopping coefficients are distinct when epsilon is nonzero. Total particle number and each pair's charge remain conserved; full flavor permutation symmetry is broken. This is not a test that removes every symmetry.

For each nonzero e, set s=sign(e), p=aL-s aR, and q=aL†+s aR†. The CAR identity is

\[
-|e|s T=-|e|I+\frac{|e|}{2}(p^\dagger p+q^\dagger q).
\]

Appending those squares and shifting b by `-sum(abs(e_i))` reproduces the perturbation exactly. The source certificate's residual is unchanged, including its norm **8.759890352758824e-7**. A nonzero residual is allowed: the certified lower endpoint is b minus that norm.

The orbit verifier now accepts an optional `direct_squares` list outside its permutation average. It applies the same rational-weight, polynomial-degree, and fixed-charge validation to both kinds of squares. No symmetry of the new Hamiltonian is assumed by this positivity argument.

## Accepted intervals

All values are in the matched model's energy units.

| epsilon | Residual-only width | Local-square width | Local-square lower | Independent upper |
|---|---:|---:|---:|---:|
| 1/1000 | 0.01199970638 | **0.00599970638** | 3.275005071503115 | 3.281004777883403 |
| 1/100 | 0.11980953318 | **0.05980953318** | 3.221005071503115 | 3.280814604687698 |

The residual-only baseline charges both hopping monomials separately. The explicit positive decomposition improves that enclosure by 0.006 and 0.06 respectively, up to the tiny existing coefficient residual. It does not recover the near-micro-unit accuracy of the symmetric source proof. This is a representation-format fix and a rigorous transfer baseline, not evidence that asymmetry has been solved computationally.

For each perturbed H, numerical diagonalization proposes an integer upper vector in ascending fixed-N bitstring order. A separate standard-library replay recomputes its exact rational Rayleigh quotient using CAR signs and the certificate's own Hamiltonian. Stored numerical energies, norms, and upper values are not trusted. The upper proposal is the only part using a many-body eigenvector; it does not supply the lower proof's directions.

Generation and both exact checks took 8.67 and 9.07 seconds in the recorded runs. No asymptotic performance claim follows from those timings.

## Reproduction and corrected diagnostics

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_symmetry_transfer --epsilon 1/1000
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_symmetry_transfer --epsilon 1/100
python3 -S -m experiments.marginal_symmetry_transfer --verify results/marginal_symmetry_transfer/hopping_1_1000/certificate.json
python3 -S -m experiments.marginal_symmetry_transfer --verify results/marginal_symmetry_transfer/hopping_1_100/certificate.json
```

Early diagnostic artifacts are preserved under `results/marginal_symmetry_transfer/invalidated_diagnostics/` and must not be used as accepted receipts. They mapped a symmetric amplitude recipe incorrectly into the Fock basis, incorrectly classified nonzero residuals as uncertified, and reported an incorrect residual scale. The corrected implementation uses a fresh integer vector and exact residual enclosure. No conclusion about a sign failure or a fundamental inability to represent the perturbation survives those corrections.
