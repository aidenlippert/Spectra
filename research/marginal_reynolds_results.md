# A small coefficient representation for the matched quartic problem

The complete quartic search on the ten-mode, five-particle matched model now has a replayed rational interval

```text
3.281005630613246 <= E0 <= 3.281006695574010
width = 1.064960764376076e-6
```

The Hamiltonian uses interaction coefficient one and hopping coefficient 1/5. These are model energy units, not Hartree molecular energies. The independent sector reference brackets the exact physical energy between 3.2810066955735238 and 3.2810066955744333.

## The representation that worked here

Work in the real Hermitian CAR polynomial algebra, quotient out the fixed-number ideal, and average over the matched model's flavor-permutation and global left/right symmetry. The numerical problem then sees the symmetry-invariant coefficient coordinates rather than every repeated monomial coefficient.

For any exported result, the exact checker reconstructs

\[
H-bI=\sum_j p_j^\dagger p_j+(\widehat N-N)X+R.
\]

On the declared particle-number sector the ideal term vanishes, each square is positive, and every normal-ordered monomial has operator norm at most one. Therefore `b - sum(abs(R_coefficients))` is a rigorous lower bound. Discovery may use floating spectral decompositions and a numerical solver; acceptance uses rational arithmetic only.

`experiments/marginal_reynolds.py` implements three reductions: charge blocks, stabilizer spectral subspaces, and full Reynolds averaging across equivalent charge blocks. The last step was decisive for this calculation.

| Formulation | Coefficient equations | PSD scalar variables | Largest PSD block | Accepted interval width |
|---|---:|---:|---:|---:|
| Complete quartic before symmetry reduction | 1,531 | 115,221 | 186 | no completed accepted raw solve |
| Stabilizer splitting | 1,531 | 18,673 | 14 | 0.0051409913 |
| Full Reynolds reduction | 39 | 2,082 | 14 | 0.0000010650 |

The final problem has 180 PSD blocks, 30 flavor-charge orbits, and 19 retained number-multiplier directions. Assembly took 23.22 seconds; the SDP solve took 1.15 seconds. These are separate timings, and the solve timing excludes assembly and certificate export/replay. The ordinary exported certificate is approximately 1.8 MB.

The solver's objective was 3.2810066958660005 and its status was `optimal`. Its smallest raw Gram eigenvalue was approximately -3.89e-12. Neither the objective nor the status supplies the rigorous bound; exact expansion charges all projection, clipping, and rounding errors. The resulting residual L1 is 1.065252754203e-6.

## What this establishes

The cubic obstruction is real, as shown by the separately checked nonrepresentable moment witness. A quartic certificate closes this instance to approximately one part in a million of the chosen energy unit. The earlier 0.00514 interval is not a lower limit on what quartic certificates can achieve.

This is not an exact proof of equality between the quartic optimum and the physical energy, or a theorem about all particle counts. The finite exact symmetry argument admits an averaged optimum; floating spectral clustering and QR rank decisions are numerical proposals and do not themselves certify losslessness.

Symmetry reduction of SOS programs is established methodology: [Gatermann and Parrilo, *Symmetry groups, semidefinite programs, and sums of squares*](https://arxiv.org/abs/math/0211450). Here it is applied to fermionic CAR coefficient maps with an independent exact checker. The implementation still explicitly enumerates a group of order 240 and expands the certificate. Its present group-enumeration strategy grows factorially with the number of matched flavors; the measured small solve is not a demonstrated large-system scaling law.

## Replay

```sh
python3 -S -m experiments.marginal_symbolic --verify results/marginal_reynolds/m10_d4_exact/certificate.json
```

The upper witness is bound to the matched Hamiltonian and independently evaluated with the rational symmetric Rayleigh routine in `experiments/marginal_collective.py`.

## Generator-based coefficient projection

A later implementation replaces explicit group enumeration in the coefficient projection with breadth-first signed orbits of adjacent flavor swaps and the global left/right swap. Direct comparison with full polynomial averaging passes. The invariant row count stays exactly 39 at 8, 10, 12, and 14 modes; projection construction took approximately 0.05, 0.25, 0.93, and 2.68 seconds respectively. These measurements concern coefficient projection only. Stabilizer construction and Gram export still enumerate the symmetry group. The uniform 39-coordinate count and 20-dimensional truncated quotient are derived in `research/marginal_invariant_dimension.md`.
