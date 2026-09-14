# Exact spectral residual witnesses

For a balanced normal-ordered k-body residual, index the real Hermitian
coefficient matrix by increasing k-tuples. Ascending annihilation order gives
the coefficient conversion (-1)^(k(k-1)/2). The dimensions are binomial(M,k):
individual orbital indices are incorrect when k exceeds one.

For each exact connected component of this coefficient matrix V, numerical
eigendecomposition proposes a rational shift ell and integer factor F/d.
The independent standard-library replay forms

    D = V - ell I - F F^T / d^2

using a common integer denominator. If g is the smallest Gershgorin endpoint
of D, then V >= (ell+g) I. The proof trusts only integer factors, never a
claimed Gram matrix or floating-point eigenvalue. It checks symmetry, factor
shape, denominator, full tuple coverage, and source certificate hash.

Take the minimum endpoint across components, including zero rows. The
exterior-power identity sum_I A_I^dagger A_I = binomial(Nhat,k) lifts this to
R_k >= binomial(N,k) min_component(ell+g) on fixed total N. Take the maximum
with the previous coefficient-L1/Gershgorin bound for each body. Handle the
constant exactly; k>N acts as zero. Sum these corrections with the original
SOS scalar b. There is no assumption on the spin sector.

The proposal and replay gate each component at dimension 600 and aggregate
cubic work at 150 million products. This is a fixed-body-rank polynomial
coefficient calculation, not enumeration of binomial(M,N) many-body states.
It does not improve the asymptotic dictionary-discovery cost. All extra
factor entries, bytes, and replay time must be added to certificate costs.

The first agent draft collapsed tuple indices, omitted the phase, and summed
component minima. Its base `wedge_spectral/receipt.json` is explicitly INVALID.
The root replacement and its H6/H8 subdirectory artifacts supersede it. An
independent read-only audit found the replacement sound within this scope.
Five exact tests cover indefinite matrices, arbitrary proposed factors,
malformed data, k=2 tuple structure, k=3 sign and particle-count lifting, and
minimum-over-components aggregation. Full molecular replay is performed with
`python -S`.
