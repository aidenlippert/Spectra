# Exact asymmetric transfer through spin-sector Schur complements

The unequal-hopping M=10,N=5 test now has a certified interval of width **2.8389575005161833e-7** at epsilon=1/1000. At epsilon=1/100 the width is **4.195706286226551e-5**. Both bounds use exact rational positivity and an independently evaluated integer-vector upper witness. The proof exploits the surviving pair-charge algebra and the matched Hamiltonian's spin decomposition. It does not solve the general marginal cone or determine the optimum of the quartic relaxation.

| Formulation | Width, epsilon=1/1000 | Width, epsilon=1/100 |
|---|---:|---:|
| Six-dimensional Schur bound | 6.764661744689823e-6 | 8.951825026491765e-4 |
| Exchange parity, three-dimensional bound | 2.550450377619148e-6 | 3.070379357854013e-4 |
| Spin-resolved nested Schur bound | **2.8389575005161833e-7** | **4.195706286226551e-5** |

The final intervals are

```
epsilon=1/1000: 3.2810044939876533 <= E0 <= 3.2810047778834033
epsilon=1/100:  3.2807726476248353 <= E0 <= 3.2808146046876976
```

The implementation is `experiments/marginal_schur_transfer.py`. Accepted certificates, search endpoints, receipts, and independent standard-library replays are under `results/marginal_schur_transfer/1_1000_resolved/` and `1_100_resolved/`. Earlier six- and three-dimensional formulations are retained in separate directories.

## The algebra that makes the bound small

Each pair has conserved charge nL+nR. In a singly occupied pair, T+T† acts as sigma_x. Empty and doubly occupied pairs have zero hopping. On every half-filled charge sector,

\[
H_0=15/4+J_z^2-(2/5)J_x.
\]

Five singly occupied pairs contain spins j=5/2,3/2,1/2. Sectors with one double and one empty pair have at most j=3/2; sectors with two of each have j=1/2. Their ground energies are bounded by the already verified small Jacobi matrices in `marginal_sector_reference.py`. Spin multiplicities do not change the block energies. Tests independently reproduce the full 252-dimensional fermionic spectrum from these blocks with multiplicities 1,24,75.

Write V=-sum_i(i-2)sigma_x(i), so H_epsilon=H0+epsilon V and ||V||<=6. For the fully symmetric spin-5/2 subspace P,

\[
PVP=0,\qquad PV^2P=(25I-4J_x^2)/2.
\]

These identities follow from centered coefficients and permutation averaging. The projected-square identity is checked against all 32 five-spin bitstrings. In the symmetric-function basis f(k), the metric is D=diag(binomial(5,k)). The matrices become symmetric only after multiplication by D; rational LDL checks this symmetry explicitly.

Global spin flip F=product(sigma_x) commutes with both H0 and V. On a spin-j component formed from five spins it acts as (-1)^(5/2-j) times Dicke reflection: each singlet contributes -1. Consequently, in the F=+ sector:

- P is the reflection-even spin-5/2 space, dimension three.
- R1 is the reflection-odd spin-3/2 space, including its multiplicities.
- R2 is the reflection-even spin-1/2 space, including its multiplicities.

The R2 energy is exactly c2=19/5. The weighted two-dimensional R1 block has matrix [[6,-3/5],[-1/5,22/5]] before multiplication by metric diag(2,6). Exact LDL proves its energy is greater than **c1=216411/50000=4.32822**.

The other pair-charge sectors and all F=- states are bounded below by 7/2 before perturbation. The verifier requires **7/2-6|epsilon|>b**, ensuring they cannot violate the claimed global lower bound.

## Why the second Schur step helps

The perturbation maps P only into R1. There is no direct P-to-R2 coupling. The verifier checks this selection rule exactly: for the symmetric bitstring columns Z, J²VZ=(15/4)VZ. It uses the identity J²=-5/4 I+sum_(i<j)Swap_ij, so the check is integer arithmetic over 32 auxiliary spin bitstrings.

For a proposed lower endpoint b, define

\[
q_2=c_2-b-6|\epsilon|,
\qquad
q_{\rm eff}=c_1-b-6|\epsilon|-36\epsilon^2/q_2.
\]

Both must be positive. Eliminating R2 first bounds the remaining R1 operator below by q_eff. Eliminating R1 then gives the sufficient condition

\[
P H_0 P-bP-\frac{\epsilon^2}{q_{\rm eff}}PV^2P\succ0.
\]

The final condition is a three-dimensional weighted rational LDL test. The term 36 comes from ||V||²<=36. No numerical ground vector is needed for this lower-bound proof. The positivity conditions are conservative; the failed bisection endpoint is a failed sufficient condition, not an upper bound on the physical ground energy or an exact cone optimum.

## Validation and practical scope

The final independent replays ran with `python -S`. They reconstruct the exact model Hamiltonian from epsilon, recompute every sector and Schur positivity gate, validate the integer upper witness, and evaluate its Rayleigh quotient using CAR signs. Stored numerical energies and search metadata are not trusted.

```sh
python -S -m experiments.marginal_schur_transfer --verify results/marginal_schur_transfer/1_1000_resolved/certificate.json
python -S -m experiments.marginal_schur_transfer --verify results/marginal_schur_transfer/1_100_resolved/certificate.json
```

Recorded generation times were about 0.15–0.19 seconds, including lower search, numerical proposal of the 252-component upper witness, and exact replay. These are individual local timings, not a controlled benchmark or a scaling result. The final positivity matrix has dimension three, but verification also uses small complementary matrices, a 32-bitstring selection-rule check, and the 252-component upper witness.

The surviving pair-charge structure is essential. General hopping between pairs, unrestricted molecular integrals, growing-system cost, and automatic discovery of the useful spin decomposition remain outside this certificate. This result supplies a strong asymmetric finite-model bound without claiming that the quartic SOS relaxation has converged.

A subsequent [general perturbation extension](marginal_general_schur_results.md) now permits the new Hamiltonian to break pair charges and exchange symmetry. It uses the same reference gap with automatically discovered coupling-response recurrences. That extension retains a specially structured reference and perturbative accuracy; it does not change the scope of the spin-specific certificates above.
