# Measured acceleration campaign

The original H6 run spent 300.95 s in stage 1 and 12.98 s in stage 2,
out of approximately 410.67 s overall. A short cProfile probe of the frozen
CPU implementation attributed 5.88 s to PSD projection (5.39 s inside eigh),
2.43 s to initial QR, and 0.80 s to sparse matrix/vector products.
Profiling overhead means those fractions are diagnostic, not benchmark timings.

Ranked hypotheses, before GPU solver changes:

1. Grouping equally sized PSD eigensolves on the A100 reduces their cost.
   Falsify by synchronized FP64 timings including host/device transfers.
2. A different CPU LAPACK eigensolver or concurrent independent blocks can
   outperform GPU transfer overhead for these small matrices.
   Compare against identical matrices and a single-thread CPU reference.
3. Keeping coefficient products and iterates resident on the GPU could remove
   hybrid transfers, but sparse solves and small-kernel launch costs could lose
   that benefit. Only pursue if timings support it.
4. The fixed 300 s first stage can do work beyond the original 1.6 mHa target.
   An earlier refinement stage may improve time to verified accuracy. This is
   a separate algorithmic experiment, never counted as hardware-only speedup.

All numerical candidates require the original rational accepting checker.
Prepared Hamiltonian maps, upper-state data, and nonsinglet proof are frozen
inputs to this component experiment; their original construction costs are
not erased from any complete-runtime comparison.

One newly created A100 SXM4 is owned by this campaign at $1.99/hour.
The two existing A10 instances are outside this campaign and remain untouched.
