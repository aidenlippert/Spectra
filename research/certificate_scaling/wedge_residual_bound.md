# Wedge-sector residual bound for exact CAR remainders

Let a balanced normal-ordered residual of body rank `k` be

`R_k = sum_(I,J, |I|=|J|=k) V[I,J] a_I^dagger a_J`,

where `I,J` are increasing mode tuples and `V` is Hermitian.  This convention
has no hidden factorial: each increasing tuple is one CAR monomial.  On the
`k`-particle wedge, `V` is the coefficient matrix in the basis
`a_I^dagger|0>`.  With canonical annihilation in ascending tuple order, the
common sign `(-1)^(k(k-1)/2)` is absorbed consistently into the definition of
the wedge basis and does not change its eigenvalue bound.

If `V >= l_k I` on `wedge^k`, then on the fixed `N`-particle sector,

`R_k >= l_k * binom(N,k) I`.

Indeed, the lifted `k`-body operator is the `k`th exterior-power second
quantization.  The identity

`sum_(|I|=k) a_I^dagger a_I = binom(Nhat,k)`

holds on every occupation determinant: exactly `binom(n,k)` increasing
`k`-subsets are occupied when `Nhat=n`.  Therefore the positive operator
`R_k-l_k sum_I a_I^dagger a_I` lifts to a positive operator on every sector,
and the displayed fixed-`N` bound follows.  Terms with `k>N` vanish.

For a cheap exact rational lower bound, use Gershgorin:

`l_k = min_I (V[I,I] - sum_(J != I) |V[I,J]|)`.

Zero rows are included and contribute `0`; omitting them would be unsound.
Every quantity is rational if the residual coefficients are rational.  The
resulting alternative lower endpoint for a certificate with scalar part `b`
and any exactly represented ideal terms is

`b + R_constant + sum_(k <= N) binom(N,k) l_k`.

Taking the maximum of this endpoint and the existing coefficient-l1 norm
endpoint is valid because both are independently certified lower bounds.

This argument uses only fixed total particle number.  It does not assume a
fixed spin or `S_z` sector.  If additional conserved occupations are used,
the matrix must be block diagonal in those labels and the bound must be
checked on every block.  For a decomposition into alpha-occupation blocks,
one may instead use the exact blockwise lower bound

`min_{N_alpha} sum_k l_k * binom(N_alpha,k) * binom(N-N_alpha, k-k_alpha)`,

with the appropriate split of the `k` occupied modes; this refinement is
valid only after proving the block decomposition and checking all integer
allocations.  It must not be inferred from an `S_z` label alone.

The bound is a residual alternative, not a claim about certificate discovery.
For fixed body rank, sparse Gershgorin rows are cheap to verify; the number of
wedge indices can still grow combinatorially with the mode count.  No
generic polynomial-size theorem follows from this observation.
