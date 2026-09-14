# V11 exact matrix arithmetic preflight

The first case is fixed in the preceding V10 note: n3 XXZ, gamma2, T1/2,
central Z and operator-norm tolerance1/1000. The current adaptive V8
fraction-free Taylor implementation and the matrix implementation must export
exactly equal polynomial coefficients and the same residual norm witness.
The original V7 checker remains unchanged. No weaker error bound, damping
change, heldout input or acquired method is present.

One initial call to each arm verifies equality and warms imports/caches.
Record these costs separately; do not use them for a speedup claim. Then run
seven pairs, shuffling arm order with seed11 within each pair. Include model
and generator creation, all construction and norm-proposal attempts, a fresh
original checker and exact endpoint evaluation. Report each phase and the
complete time. No norm witness is rebuilt unnecessarily outside the adaptive
method. JSON serialization is outside the calculation timer.

Matrix arithmetic retains real and imaginary integer arrays. Before int64
kernels, forward bounds must cover temporary values, not just final cancelled
results. A costed Python-integer object-array fallback handles wider values.
The matrix route is bounded to n<=6 and still exports at most512 Pauli terms
per coefficient; the fresh checker enforces its unchanged cache/live caps.

A loss closes this candidate as currently implemented; it does not establish
optimality of the sparse implementation. A win would only justify expansion
to development cases and a stronger cost study. Supplied matrix/FWHT identities
cannot be labelled m1, and no learner or reserved evaluation follows from an
isolated faster kernel.
