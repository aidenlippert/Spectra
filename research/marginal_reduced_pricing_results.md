# Stabilizer pricing without large moment eigensolves

The matched ten-mode search now has an optional `--pricing stabilizer` path. It replaces 30 full candidate eigensolves, of maximum dimension 186, with 180 subblock eigensolves of maximum dimension 14. Subblocks of equal dimension are contracted and diagonalized in batches.

This is a reduction of repeated pricing work. Preparation still constructs the full quartic word dictionaries and coefficient maps, and enumerates stabilizer permutations. Those scaling problems remain.

## Why all negative directions remain detectable

For a charge dictionary, let `G` map Gram coefficients to the Reynolds-invariant coefficient rows. The priced moment matrix is the symmetric part of `reshape(G.T @ dual)`. It commutes with the signed permutation representation of the dictionary's stabilizer, because the coefficient functional is invariant.

`split_blocks` diagonalizes a real symmetric linear combination of stabilizer representation matrices. The moment matrix commutes with this combination and therefore preserves each of its eigenspaces. In exact arithmetic its spectrum is the union of the compressed spectra. Degeneracy can leave a larger subspace; it does not require dropping any directions. All subspaces are retained.

The implementation uses floating eigenspace bases and numerical clustering. Tests therefore reconstruct the full moment matrix and compare its complete spectrum on six-mode dictionaries. The ten-mode benchmark compares negative dictionary sets, minimum eigenvalues, and lifted columns on 20 fixed duals. These are finite numerical checks, not exact certification of the floating decomposition. Final energy certificates still pass through the rational algebra checker independently.

## Measured pricing cost

Five alternating timing batches, each using the same 20 dual vectors, gave:

| Path | Median seconds per pricing call |
|---|---:|
| Full moment matrices | 0.00469335 |
| Batched stabilizer subblocks | 0.00094501 |

This is approximately **4.97 times faster for the measured pricing operation**, including map contraction, eigenvectors, and selected column construction. It is not a factor of 4.97 for the entire search. Preparing the reduced path took 3.82 seconds in this run.

Across the benchmark duals, corresponding minimum eigenvalues differed by at most 4.89e-15 and lifted columns by at most 8.33e-16. The initial unbatched implementation was slightly slower than full pricing; its receipt is preserved to distinguish the effect of batching from matrix dimensions alone.

Eigenvectors in a degenerate minimum eigenspace are not unique. The two methods can choose different valid cuts, so equal minimum eigenvalues do not imply identical adaptive trajectories or identical accuracy at a fixed cut budget.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_pricing_benchmark
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_energy_selector --strategy negative --iterations 117 --max-cuts 117 --pricing stabilizer
```

The benchmark receipt is `results/marginal_reduced_pricing/m10_benchmark.json`. The original full pricing path remains available as the default for comparison.

## Complete discovery run

The stabilizer pricing run added 117 scalar cuts and independently checked the interval `3.281004872488434 <= E0 <= 3.281006695574010`, width **1.8230855759971247e-6**. Its numerical objective was 3.2810049683237636 and its exact residual norm 9.56e-8. Total observed time was 42.22 seconds, including export and rational checks. This run does not establish an end-to-end speedup: the earlier full-pricing run at the same cut budget produced the tighter width 8.157e-7, and the timings were not controlled end-to-end comparisons. The retained default remains full pricing.
