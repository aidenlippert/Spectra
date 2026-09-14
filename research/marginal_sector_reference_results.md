# Locating the ten-mode certificate gap

The matched ten-mode, five-particle Hamiltonian at t=1/5 has ground energy in the exact interval

\[
\frac{1803752506297}{549755813888}
< E_0 \le
\frac{3607505012595}{1099511627776}.
\]

The width is 1/1099511627776, approximately 9.095e-13. Numerically E0 is approximately 3.281006695574. The previous rational variational upper witness exceeds the true ground energy by at most 4.865e-13. The earlier gap near 0.02385 therefore was not caused by an inadequate upper state.

## Reduction covering all particle sectors

Let M=2m and N=m, with matched pairs (i,i+m). Each pair occupation F_i=n_i+n_{i+m} is conserved. If d pairs are doubly occupied, exactly d must be empty, leaving s=m-2d singly occupied pairs. All choices of which pairs are empty or double have identical spectra for this matched, equal-weight model.

Order fermion modes pair by pair. A matched hopping term then acts between adjacent modes, with matrix element -t on its singly occupied pair. It does not acquire a spectator-dependent fermion sign. The active part of a fixed sector is an s-bit hypercube. A vertex with k left occupations has potential

\[
V_d(k)=\binom{d+k}{2}+\binom{d+s-k}{2}.
\]

For t>0 and s>0, this real matrix has negative hopping entries on a connected graph. Apply Perron–Frobenius to a sufficiently shifted negative of the Hamiltonian: its ground vector is unique and strictly positive. Permuting active pairs commutes with the matrix, so uniqueness and positivity force that vector to be permutation invariant. The ground state therefore lies in the normalized Hamming-weight basis, whose Jacobi matrix is

\[
(J_d)_{kk}=V_d(k),\qquad
(J_d)_{k,k+1}=-|t|\sqrt{(k+1)(s-k)}.
\]

A bipartite sign transformation covers negative t. At t=0 the diagonal minimum is represented in the same symmetric space; s=0 is a scalar block. Thus, including all d rather than assuming the singly occupied sector wins,

\[
E_0=\min_{0\le d\le\lfloor m/2\rfloor}\lambda_{\min}(J_d).
\]

This exact reduction relies on matched equal hopping and the two uniform density-interaction groups. It is not valid for the asymmetric extra-edge control without a different derivation.

## Rational spectral check

Square roots need not be approximated. For a rational trial bound b, compute leading determinants of J_d-bI:

\[
D_0=1,\quad D_1=V_d(0)-b,
\]
\[
D_{k+1}=(V_d(k)-b)D_k-t^2 k(s-k+1)D_{k-1}.
\]

Sylvester's criterion says all leading determinants are strictly positive exactly when b lies strictly below the smallest eigenvalue. Bisection therefore gives rational endpoints. A zero determinant is treated as failure of strict positive definiteness, which safely places the trial point at or above the smallest eigenvalue.

The initial lower endpoint is min(V_d)-|t|s-1, using the hypercube hopping norm |t|s. The initial upper endpoint min(V_d) is the energy of a basis state. The implementation checks both endpoint predicates before bisection.

For ten modes the entire ground-energy comparison uses matrices of dimensions 6, 4, and 2:

| Doubly occupied pairs d | Active pairs s | Ground energy, approximately |
|---:|---:|---:|
| 0 | 5 | 3.281006695574 |
| 1 | 3 | 3.551000400320 |
| 2 | 1 | 3.800000000000 |

Their rational brackets are disjoint. The d=0 sector is rigorously the lowest for this instance.

## Completing the pure-triple dictionary

The prior full mixed-cubic dictionary included every mixed a_k† a_j a_i word but only within-half pure triples. A new experiment replaces those pure triples with every annihilation triple and its adjoint, retaining the full mixed dictionary and exact charge splitting.

| Cubic dictionary | Numerical proposed b | Exactly verified lower bound |
|---|---:|---:|
| Edge subset, prior run | 3.257178392900 | 3.25716083471628 |
| Full mixed, local pure triples, prior run | 3.257178366022 | 3.25715679175004 |
| Full mixed, all pure triples, new run | 3.257178393168 | 3.25716059376644 |

The new residual allowance is 0.00001779940156. Its certified interval width remains 0.02384610180757. Completing the pure triples did not close the gap in this numerical run.

The reference proves that the existing certified lower bounds fall short of the physical ground energy. Agreement of three numerical proposals supports a limitation of the current relaxation, but does not prove its mathematical optimum. No exactly feasible dual moment witness or other upper bound on the relaxation optimum has been certified. Numerical optimization limitations remain a separate unresolved possibility. In particular, these results do not establish that every cubic SOS certificate must fail.

The next decisive step is to certify the relaxation's dual optimum, or to construct a stronger certificate that closes the gap. The exact spectral reference now allows either route to be evaluated without uncertainty from the upper-state ansatz.

## Implementation and validation

`experiments/marginal_sector_reference.py` uses only standard-library rational arithmetic and does not enumerate Fock states. It can bind the reference calculation to an existing certificate by checking equality with the exact matched Hamiltonian and the required particle count.

`experiments/marginal_cubic_probe.py` discovers and exports the expanded cubic certificate using the existing numerical solver and exact checker. It uses no many-body diagonalization.

Four new tests check the spectral predicate, argument and Hamiltonian rejection, complete pure-triple coverage, and agreement with an independent full fermion-matrix calculation for M=4,6,8,10 at t=0,1/5,-1/3. This independent matrix enumeration occurs only in tests. All 37 tests across the marginal experiment modules passed.

The new certificate and reference replay successfully in `python3 -S`, with NumPy, SciPy, and CVXPY confirmed absent from the process. Results and source hashes are in `results/marginal_sector_reference/`.

```sh
python3 -S -m experiments.marginal_sector_reference --certificate results/marginal_compression/m10_t1_5_full_split.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_cubic_probe
python3 -S -m experiments.marginal_collective results/marginal_sector_reference/complete_cubic.json
```
