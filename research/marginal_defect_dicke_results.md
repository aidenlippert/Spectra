# Implicit local-defect algebra for exact Schur certificates

The lower-bound calculation can now avoid explicit configuration vectors on four existing certificates. A new local-defect/Dicke representation reproduces both six-dimensional reference-resolvent proofs and both first-round enlarged proofs, including the 36-dimensional mixed-interaction construction. Every bound is accepted by rational positivity checks.

This removes configuration enumeration from these lower-bound constructions. The four replays in this report retain explicit upper witnesses. A subsequent [complete implicit construction](marginal_implicit_certificate_results.md) now ports targeted enrichment and variational upper witnesses, with both tested intervals below 1e-7.

## The representation

Arrange m matched pairs in the Fock gauge L0,R0,L1,R1,... . An atom fixes an occupation o_i in {empty,L,R,LR} on a distinguished set S of pairs. Every other pair is singly occupied, and the atom sums all assignments with exactly k right occupations among those spectators. Denote it by |S,o;k>.

The original unnormalized symmetric vector Z_k is the atom with S empty. A general represented state is a rational linear combination of atoms. Different atom combinations can represent the same physical vector: this is an overlapping generating family, not an independent basis.

### Exact overlaps

For two atoms, active occupations must agree on their shared distinguished pairs. An occupation distinguished in only one atom must be singly occupied; otherwise the overlap is zero. Subtract the right occupations fixed by the other atom from each spectator right count. If the two remaining counts disagree, the overlap is zero. Otherwise, for their common remaining count K,

\[
\langle S,o;k\mid T,p;\ell\rangle
=\binom{m-|S\cup T|}{K}.
\]

Out-of-range counts give zero. This count requires no spectator configuration enumeration. It also detects exact physical null relations that would be missed by treating atom coefficients as independent coordinates.

### Fermionic actions

Canonical input modes L0,...,L(m-1),R0,...,R(m-1) are mapped to the pair-contiguous gauge. Applying a CAR operator to an undistinguished pair branches only over its two possible single occupations and makes that pair distinguished. Each right branch reduces the remaining spectator right count by one.

The sign exponent for an operator on pair i is

\[
i+\sum_{j\in S,\,j<i}(\operatorname{popcount}(o_j)-1)
+\mathbf1_{\text{right mode}}\mathbf1_{\text{left occupied}}.
\]

The baseline i counts one electron in each earlier spectator pair; distinguished occupations correct that count. Tests independently expand atoms into canonically ordered Fock vectors and compare single operators, cross-pair products through degree eight, and CAR anticommutators.

### Reference action

The matched reference H0 never enlarges S. Its interaction energy is the scalar binom(NL,2)+binom(NR,2), determined by the active occupations and k. Hopping flips active L/R occupations. On n spectators it acts as

\[
\sum_i a_{R_i}^\dagger a_{L_i}|D_k^n\rangle
=(k+1)|D_{k+1}^n\rangle,
\]

with reverse coefficient n-k+1. Thus arbitrarily many reference actions preserve the distinguished-pair support of every atom.

## Exact recurrence and retained-space construction

The coupling W=(I-P)(H-H0)Z is built and projected using these overlaps. An exact scalar Lanczos recurrence in the Hilbert space of block columns discovers an annihilating polynomial. Termination requires the residual's exact squared norm to be zero. The projected Hamiltonian, coupling Gram matrix, annihilator, and every response moment agree exactly with the explicit calculation on both fixtures.

For enlargement, a Gram-LDL procedure selects independent H0-Krylov vectors. It admits a new vector only when its exact residual norm is positive, and rejects dependent candidates only at zero exact norm. Every admitted vector's H0 image is processed, establishing closure without a configuration basis.

For the retained columns U, compute G=U^TU, A=U^THU, and B=(HU)^T(HU). The leakage Gram matrix follows from

\[
L^TL=B-A^TG^{-1}A.
\]

The existing exact Schur test then applies in this independently selected nonorthogonal basis. The checker verifies that leakage from the initial reference columns vanishes exactly.

## Replayed certificates

| Certificate | Matrix dimension | Existing certified width | Implicit representation |
|---|---:|---:|---|
| Cycle, strength 1/1000, reference resolvent | 6 | 5.603353530955592e-7 | 80 coupling atoms, at most 2 distinguished pairs |
| Mixed, strength 1/1000, reference resolvent | 6 | 5.996685860823592e-7 | 100 coupling atoms, at most 3 distinguished pairs |
| Cycle, strength 1/100, first enlargement | 14 | 2.116834216764047e-5 | 170 retained atoms; 1,181 action atoms |
| Mixed, strength 1/100, first enlargement | 36 | 2.1972817593562723e-5 | 1,912 retained atoms; 12,122 action atoms |

Atom totals count terms across all columns, including repetitions. The action columns touch at most four distinguished pairs per atom for the cycle and five for this mixed fixture. The first enlarged implicit replays took 0.78 seconds and 100.07 seconds in individual runs. A profile of the cycle replay attributed 1.058 of 1.207 instrumented seconds to inner products. The implementation now clears denominators per vector, accumulates integer overlap sums, and computes only one triangle of symmetric Gram matrices. Subsequent standard-library runs took 0.49 seconds and 54.13 seconds; every exact positivity pivot and interval endpoint matched the pre-optimization run. These individual runs were not an isolated performance benchmark. The mixed implicit construction remains slower than its earlier explicit counterpart: eliminating configuration enumeration does not yet establish an overall speedup.

An additional diagnostic applies one fixed two-pair Hermitian hopping operator and seven H0 powers at 10,20,40,100,200 modes. At 200 modes it uses at most 30 atoms, each with two distinguished pairs. Exact squared norms are stored. This is a finite local-response calculation, not a 200-mode energy certificate or a test of dense arbitrary perturbations.

## What is proved about size

After r applications of a degree-at-most-four perturbation, each atom has at most 4r distinguished pair labels. Intervening H0 actions do not change that bound. Consequently the total possible atom dictionary is bounded by

\[
\boxed{\sum_{q=0}^{\min(4r,m)}\binom mq4^q(m-q+1).}
\]

For fixed r this is O_r(m^(4r+1)). The bound is for the complete generating dictionary, not a rank formula or a claim that the displayed exponent is optimal. At N=m, a valid active pattern must have exactly q particles, equivalently equal numbers of empty and doubly occupied active pairs. The spectator right count remains free; particle number does not fix it.

This establishes a polynomial atom-count bound at fixed perturbation order around this reference. It does not bound how large r must become for a requested accuracy, furnish an efficient implicit representation for a general reference, or establish favorable arithmetic cost. Overlap assembly, redundancy, and rational coefficient growth still matter.

## Remaining work and validation

Targeted enrichment and implicit upper witnesses are now implemented in the subsequent construction, and charge-sector grouping reduces overlap assembly. Energy certificates at larger sizes remain a next step. General reference discovery and accuracy control when the reference gap is consumed remain separate mathematical problems.

Implementation is `experiments/marginal_defect_dicke.py`. The complete marginal suite passed **129 tests in 33.295 seconds**, including nine new tests. These cover overlapping atoms and null relations, fermionic signs, fast H0 action, exact explicit-versus-implicit coupling and moment agreement, closure rank, certificate replay, malformed upper witnesses, and false lower bounds. The 100-pair test checks binomial norms and local actions without expanding configurations.

```sh
python3 -S -m experiments.marginal_defect_dicke --verify results/marginal_general_schur/mixed_1_1000_paired_resolvent/certificate.json
python3 -S -m experiments.marginal_defect_dicke --verify results/marginal_enlarged_schur/mixed_1_100/certificate.json
```

Receipts and the local-response diagnostic are under `results/marginal_defect_dicke/`. All four implicit paths passed independent standard-library replay. Both enlarged paths also passed standard-library replay after integer accumulation was introduced, with every rational pivot matching the earlier output. Original certificates are preserved. These results extend [the enlarged Schur construction](marginal_enlarged_schur_results.md); they do not establish the general physical-marginal cone's boundary.
