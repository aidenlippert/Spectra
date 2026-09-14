# Compressing certificates by conserved charges

The eight-mode mixed-cubic Gram blocks can be split from size 232 to size 28 without changing the mathematical relaxation. This follows from conserved component occupations of the actual hopping graph. A further heuristic restriction to hopping-edge words retains tight bounds at three tested couplings. At ten modes, both the restricted and full dictionaries leave a gap near 0.02385. The latter result is an unresolved limitation, not a successful transfer of the eight-mode accuracy.

## Exact reduction

For each connected component C of the one-body hopping graph, define F_C=sum_{i in C} n_i. The implementation checks that every nonzero Hamiltonian monomial has zero charge under every F_C. Density interactions satisfy this condition. Interactions that violate the proposed symmetry are rejected.

Every ladder word w has a charge vector chi(w). Under U(theta)=exp(i sum_C theta_C F_C), it acquires the phase exp(i chi(w).theta). Averaging a feasible polynomial SOS identity over these phases removes precisely the Gram entries between unequal charges:

\[
Q_{uv}\longmapsto\mathbf{1}_{\chi(u)=\chi(v)}Q_{uv}.
\]

Each retained block is a principal submatrix and remains positive semidefinite. The Hamiltonian and scalar bound are invariant. The full Hermitian two-body number multiplier space is closed under this averaging, so the multiplier can be averaged too. Hence the full charge-split dictionary has the same exact feasible energy bounds as the original dictionary. This argument applies to the chosen relaxation, not the complete marginal cone.

Because the dictionary consists of individual ladder monomials, the operation simply zeros entries of a real symmetric Gram matrix. No pairing of positive and negative charge sectors or complex matrix representation is required. Every original word is retained.

For matched hopping, the conserved components are individual matched pairs. In the asymmetric six-mode control the extra edge merges two pairs, giving components {0,1,3,4} and {2,5}. The code discovers this change from the Hamiltonian rather than incorrectly imposing the matched symmetry.

## Heuristic reduction

All experiments retain the quadratic and local pure-triple dictionaries. Within the mixed-cubic blocks, three choices are compared:

- **Full:** all a_k† a_j a_i, with i<j, plus linear annihilators, and the adjoint block.
- **Density:** retain only k=i or k=j, giving density-dressed annihilators.
- **Edge:** retain the density words and words where k is connected by an actual hopping edge to i or j.

These subsets are selected before solving, without using the full certificate's factors. The subset restrictions are heuristic and need not preserve the relaxation. Charge splitting itself remains lossless for each subset because every retained word is an eigenoperator of the symmetry.

For a matched graph with M=2m and nonzero hopping, the edge words have the flavor charge of a single annihilator. Each such mixed block has 2[1+1+4(m-1)]=8m-4=4M-4 words: two choices of annihilator within a flavor, each multiplied by the identity, its partner density, or one of four one-body words on another flavor. This explains the observed largest block sizes 28 and 36. The unsplit full mixed block has M+M*binom(M,2) words. Other coefficient and multiplier costs still grow with mode count.

## Exact replay results

Energy units and Hamiltonian are defined in `marginal_coefficient_results.md`. All matched runs use N=M/2. Each table width is the difference between an exact rational lower certificate and an exact rational variational upper witness; displays below are rounded.

| Modes | t | Dictionary | Largest Gram block | Certified lower | Interval width |
|---|---:|---|---:|---:|---:|
| 8 | 0.2 | Full, unsplit (prior run) | 232 | 1.63782598541069 | 0.00003540930236 |
| 8 | 0.2 | Full, charge split | 28 | 1.63785338803138 | 0.00000800668167 |
| 8 | 0.2 | Density, charge split | 16 | 1.61050626448845 | 0.02735513022460 |
| 8 | 0.2 | Edge, charge split | 28 | 1.63784728750356 | 0.00001410720949 |
| 8 | 0.1 | Edge, charge split | 28 | 1.89084536733790 | 0.00000563831415 |
| 8 | 0.5 | Edge, charge split | 28 | 0.62770654113342 | 0.00001213559757 |
| 10 | 0.2 | Edge, charge split | 36 | 3.25716083471628 | 0.02384586085773 |
| 10 | 0.2 | Full, charge split | 36 | 3.25715679175004 | 0.02384990382397 |

At eight modes, charge splitting reduces independent symmetric Gram variables from 57,092 to 4,616. Edge pruning reduces this to 3,680. These are structural counts, not measured runtime speedups. The full split certificate contains 6,365 nonzero factor entries versus 32,556 in the prior unsplit export. Differences between the achieved full bounds arise from numerical solving and rational export; the exact relaxations are equivalent by the averaging argument.

The asymmetric six-mode full split run gives lower bound 0.54087009148456 with largest block 44, compared with 96 before splitting. This run supplies no new asymmetric upper witness.

At ten modes the full and edge numerical proposals are approximately 3.25717837 and 3.25717839. Their certified gaps are much larger than their residual allowances (about 2.16e-5 and 1.76e-5). Restoring the full dictionary therefore did not close the observed gap. Solver proposals alone do not certify the optimum of the relaxation. We have also not proved the symmetric upper ansatz is the exact ground state. These data do not yet distinguish relaxation looseness from upper-ansatz looseness or numerical optimization limitations.

## Implication for the next attack

The useful coordinates here are operator words classified by conserved charges, with density and hopping operators multiplying annihilators. Symmetry removes forbidden cross terms exactly; the edge subset removes many additional words empirically without losing observed eight-mode accuracy.

The next discriminating experiment should target the ten-mode gap: establish a sharper spectral reference using its conserved sectors, then test an enlarged certificate dictionary if the gap belongs to the relaxation. Increasing numerical precision alone is unlikely to remove the observed 0.02385 width given the much smaller residual allowances, although numerical relaxation optimality is not certified here.

This construction relies on the disconnected hopping graph and density interactions. A connected generic hopping graph may leave only total particle number, eliminating most of this compression. The current result does not establish a transferable compact representation for general molecules or materials.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_compression --modes 8 --subset full --split
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_compression --modes 10 --subset edge --split
python3 -S -m experiments.marginal_collective results/marginal_compression/m10_t1_5_edge_split.json
```

Eight new lower certificates and seven matched upper witnesses replayed successfully with standard-library-only Python. All 33 tests passed, including six new checks for component selection, rejection of symmetry-breaking interactions, preservation of dictionary words, exact equality of a projected square and its charge blocks, adjoint closure of pruned dictionaries, and a saved tight bound.

New implementation: `experiments/marginal_compression.py`. Results and provenance: `results/marginal_compression/replay.json` and `results/marginal_compression/manifest.json`. Historical certificates and implementation files were preserved.
