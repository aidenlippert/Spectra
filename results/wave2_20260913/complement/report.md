# Wave2 selected-subspace Schur complement

For an orthogonal selected subspace `P` with block matrix
`H=[A B; B* C]`, a valid lower endpoint is
`lambda_min(A)-||B||^2/(mu-lambda_min(A))` when an independently certified
`C>=mu I` satisfies `mu>lambda_min(A)`. Here `lambda_min(A)` and `mu` use
exact rational Gershgorin bounds and `||B||^2` uses the exact Frobenius square,
which upper bounds the spectral coupling norm. If `mu` fails the strict gate,
the code refuses.

The original rational H4 STO-3G fixture was reconstructed in the full fixed-N
sector (70 determinants). A four-determinant selected span `[15,30,45,51]`
gave `a=-1136573199/320000000` (about -3.55179), while full-complement
Gershgorin gave `mu=-4144209931929/1000000000000` (about -4.14421), hence
`mu<=a` and refusal. The complement required 4,356 entries. The diagonal-only
value is retained as a negative-control diagnostic only; it is not a bound.

The H6 fixture (924 determinants) with selected span `[63,95,111,119]` gave
`a=-612577410647/100000000000` (about -6.12577), `mu=-3863070939667/500000000000`
(about -7.72614), and refusal. Full complement storage/work was 846,400
entries. This is explicitly an oracle-scale diagnostic, not a scalable claim.

The selected states are determinant basis states from the Hamiltonian sector,
not exact eigenvectors. Exact CAR reconstruction is performed by replaying each
fixture monomial. A zero-residual excited-state refusal remains covered by the
Wave1 tests. No local parent or analytic H6 complement certificate was found in
this bounded run; that is the obstruction the route must solve.

The reproducible output is `run.json`; implementation is
`research/wave2_20260913/complement/selected_schur.py`.
