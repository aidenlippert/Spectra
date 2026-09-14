# Recognizing and certifying a restricted dressed quadratic family

This is a working structural control, not a solution of the molecular scaling
problem. It recognizes an unknown matching transformation from H, constructs
small certificates directly, and verifies an upper bound without an FCI vector.
The frozen H4/H6/H8/H10 molecular fixtures all refuse the structural condition
because of unsupported quartic interactions. The main goal remains active.

## Exact recognition statement

Let M >= 4 and let H be a real rational, canonical, number-conserving Hermitian
CAR polynomial of degree at most four. Form q from its constant and quadratic
terms, and let G be the graph of its nonzero off-diagonal hopping coefficients.
Require G minus any single vertex to be connected. A separate validated pure
quadratic path needs no connectivity assumption.

For a matching S, define U as the product of CZ gates exp(i*pi*n_i*n_j) on its
edges. A ladder operator on i transforms to itself times (1-2*n_j) if j is its
mate. For a hopping i to j, the external controls are

    D_ij = (S_i symmetric_difference S_j) minus {i,j}.

The word "external" is essential. Complementing S does not preserve the full
symmetric difference, but does preserve these external controls. Each observed
quartic controlled hopping must have exactly the coefficient -2*h_ij times
its canonical CAR sign. Density-density and pair-hopping quartics refuse.

Write t_ijk for the observed external control bit, with absent controls zero.
For each deleted vertex k, solve on G-k

    S_ik xor S_jk = t_ijk.

A BFS either refuses an inconsistent cycle/disconnected graph or determines
one representative r_ik per column, up to a bit g_k. Symmetry imposes

    g_i xor g_k = r_ik xor r_ki  for EVERY pair i != k.

Fix g_0=0, derive all other offsets, and check all pairs. There are at most
two symmetric candidates, related by off-diagonal complementation. Check the
matching degree bound and regenerate the entire original H exactly before
accepting either. This is a sufficient recognition procedure, not a search
over exponentially many matchings. Its graph work is O(M*|E|+M^2), with O(M^2)
bit storage, plus input validation and constant-degree polynomial regeneration.
It is deliberately incomplete outside the connectivity/degree hypotheses.

The complement ambiguity has a physical explanation: the product of CZ over
all pairs is (-1)^[N(N-1)/2] in the fixed-N sector. It commutes with every
number-conserving Hamiltonian.

## Compact lower and upper witnesses

Suppose recognition accepts H=U q U-dagger, with
q=c0+a-dagger*h*a and real symmetric h. Choose a scalar mu and decompose
h-mu*I=P-Q into positive and negative spectral parts. Linear annihilation
factors for P give a-dagger*P*a. Linear creation factors for Q give
trace(Q)-a-dagger*Q*a. Thus

    q = c0+mu*N-trace(Q) + SOS + mu*(number-N).

Conjugation transports each linear factor to a polynomial with at most 2M
monomials of degree at most three. There are at most M factors. No full cubic
word dictionary or semidefinite discovery problem is constructed.

The implementation rounds the factor rows to integers over a common D and
uses the EXACT rounded hole-row norm for trace(Q). The ideal multiplier is
mu times the IDENTITY. Every bound is then recomputed by the existing exact
CAR verifier, which subtracts the full coefficient-l1 residual, including
possible degree-six rounding residuals. A numerical eigenvalue is never the
accepting certificate.

For an independent upper witness, round M-by-N occupied orbital columns C.
Exact rational elimination checks rank and computes

    gamma = C*(C-transpose*C)^-1*C-transpose,
    upper = c0 + trace(h*gamma).

This is the energy of U applied to the normalized Slater state represented by
C. Exact regeneration of H validates transport of this upper bound. Computing
the trace via C-transpose*h*C uses O(M^2*N+M*N^2+N^3) rational operations.
No occupation-sector enumeration or many-body vector enters either witness.

The witness has O(M^2) coefficients and lower replay uses O(M^3) factor
products of bounded-degree words. Bit complexity must still be counted:
integer factor precision and the rational inverse grow with requested accuracy.
Compactness alone does not give fixed-precision accuracy at arbitrary M.

Existence of arbitrarily accurate rational witnesses follows by approximating
the spectral factors and occupied subspace. For example, with spectral-row
norms bounded by sqrt(L), rounding each row entry within 1/(2D) gives a bound
on the pre-dressing quadratic coefficient error of

    M^(5/2)*sqrt(L)/D + M^3/(4D^2).

Matching dressing amplifies its coefficient-l1 norm by at most nine. Rounding
a chemical potential in the occupied eigenvalue interval adds at most M/(2D)
to the optimal dual-energy deficit. If C approximates orthonormal occupied
columns with operator error delta<1/2, its projector error is at most
delta/(1-delta), and its energy excess is bounded by
2*N*||h||*delta/(1-delta). These are approximation arguments, not a validation
of the floating eigensolver's error. The implemented finite-precision proposer
has no universal accuracy/termination guarantee; exact replay decides each run.

## Evidence and scope

Artifacts are in results/certificate_scaling/dressed_quadratic/. The two-hub
family has paired vertices, two unmatched hubs coupled to every paired vertex,
onsite energies, and matching-pair hops. It has O(M) H terms, passes the graph
hypothesis, and has nonzero correlated-hopping quartics. Orbital labels were
randomly permuted before recognition; the original matching was not supplied
to discovery. The sweep uses 4,8,16,32,64 modes at half filling.

All five intervals pass 0.0016 Ha under a separate python -S replay. At M=64,
the width is 7.05541936e-6 Ha, with 64 factors and a 293111-byte full witness.
Proposal plus its embedded replay plus an additional in-process replay took
9.93 seconds; separate standard-library replay took 5.22 seconds. These are
single-run timings, not a fitted asymptotic theorem. Detailed stage times,
bytes, factor nonzeros, exact fractions and hashes are saved per size.

Tests independently apply ladder words to occupation masks, checking the CZ
matrix identity on the entire 4-mode Fock space and two-sided bounds against
fixed-N spectra for N=0,1,2,3,4. Corrupted H, invalid matching, noncanonical
input, malformed rational data, disconnected graph and rank-deficient upper
witness cases are refused. Altered integer factors remain sound only after
paying their recomputed residual.

This family is hidden free physics. It does not establish chemistry-relevant
compression, locality of generic dual certificates, or a small-certificate
theorem for Coulomb Hamiltonians. Primary precedents include Fendley's
[Free fermions in disguise](https://arxiv.org/abs/1901.08078), which solves a
quartic chain using nonlinear free-fermion operators, and Jones and Linden's
[Integrable spin chains and the Clifford group](https://arxiv.org/abs/2107.02184),
which constructs integrable families through Clifford transformations. The
recognition procedure above is our restricted derivation; no novelty claim
about the model family or free-fermion transport is made.
