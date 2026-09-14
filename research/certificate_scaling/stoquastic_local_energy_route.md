# Stoquastic local-energy certificates

Let H be a real symmetric matrix in a fixed configuration basis, with
H[x,y] <= 0 for x != y. For every strictly positive vector psi, the
Collatz-Wielandt inequality gives

    E0(H) >= min_x (H psi)[x] / psi[x].

Writing psi_theta(x) = exp(f_theta(x)) and
f_theta(x) = sum_a theta_a phi_a(x), the local energy is

    L_x(theta) = H[x,x] + sum_{y != x} H[x,y]
                 exp(f_theta(y) - f_theta(x)).

Every summand is a nonpositive constant times an exponential of an affine
function, so L_x is concave in theta. Therefore the hypograph constraints
    b <= L_x(theta) for every x
form a convex feasible set. Maximizing b over a bounded parameter box is a
convex optimization problem. If the features and moves have bounded scope and
the factor graph has bounded treewidth, dynamic programming evaluates
min_x L_x(theta) as a separation oracle. Its cost is exponential in treewidth
and local alphabet size, rather than total system size.

This is a restricted structural theorem, not a general fermionic certificate.
The Hamiltonian must be stoquastic after a fixed basis or gauge transform, and
the positive-factor family must contain or approximate the Perron ground
state. Fermionic hopping signs, flux, frustration, exchange, and particle
number constraints can violate the premise. Local terms can also produce a
factor graph with large treewidth.

## Three-state noncommuting example

Take

    H = [[0, -1, 0], [-1, 0, -1], [0, -1, 1/2]].

The hopping terms do not commute with the diagonal potential. With
f_theta=(0, theta, 2 theta), psi=(1, exp(theta), exp(2 theta)), the local
energies are

    L0 = -exp(theta)
    L1 = -exp(-theta) - exp(theta)
    L2 = 1/2 - exp(-theta).

The best bound in this one-parameter family is max_theta min(L0,L1,L2),
which can be found by checking pairwise intersections. The unrestricted exact
ground energy is the smallest eigenvalue of H. Equality requires the logarithm
of the Perron vector to be linear in this chosen feature; a richer bounded
local family preserves the same hypograph construction.

The practical target is therefore precise: identify molecular Hamiltonians
that are stoquastic in a useful basis and whose Perron log-amplitude has a
bounded-width factor representation. That condition would provide compact,
efficiently discoverable lower bounds. It does not imply chemical accuracy for
generic molecules, superconductors, or FeMoco.

## Positive-operator bridge

Put t_xy = -H[x,y] >= 0. For any positive psi, each off-diagonal edge has
the exact positive semidefinite identity

    t_xy [[psi[y]/psi[x], -1], [-1, psi[x]/psi[y]]].

Embedded in the x,y coordinates, this block is PSD and its diagonal
contribution is exactly the corresponding part of H-b when the diagonal
remainder is L_x(theta)-b. Thus a positive local-energy proof can be converted
to an operator certificate without forming a global Gram matrix.

When a move and the ratio psi[y]/psi[x] depend on only a bounded neighborhood,
all flips of the same local type can be grouped into conditional 2-by-2 PSD
blocks. The number of block templates is then controlled by local scope and
alphabet size. A tree-decomposition dynamic program can express
L(x)-b as nonnegative local reparameterized functions plus telescoping
separator messages. Each nonnegative diagonal local function is itself a PSD
diagonal block, while the conditional flip terms use the 2-by-2 blocks above.
The running-intersection property is the condition that makes the message
terms cancel when the blocks are summed globally.

The quantifiers matter. This construction assumes a Hermitian, detailed-
balance-compatible rate representation, a fixed stoquastic gauge, bounded
scope for both rates and amplitude ratios, and a valid tree decomposition of
that scope graph. A global fixed-particle-number constraint is an additional
factor and may require separator tables whose state count grows with the
constraint range; it cannot silently be treated as local.

For rational positive amplitude ratios, every 2-by-2 block is rational PSD and
can be replayed by exact rational LDL or a rank-one factor. This avoids
encoding transcendental exponentials in the certificate: use rational
positive-factor weights directly, then bound any approximation error in the
diagonal remainder.

As a positive control, take a periodic one-dimensional spin chain with
psi(x)=2^{-w(x)}, where w is the number of domain walls, and local Hamiltonian
terms consisting of a spin flip plus the diagonal ratio contribution

    h_i(x) = -X_i + 2^{w(x)-w(x with spin i flipped)} on the diagonal.

The ratio depends only on the two neighboring spins, so each h_i is a finite
set of neighbor-conditioned 2-by-2 PSD blocks. The terms generally do not
commute, while the chosen positive state has a known zero local-energy
certificate by construction. This is a controlled validation family for the
bridge; it is not evidence that generic fermionic molecular Hamiltonians meet
the same locality and sign conditions.

## Discovery conditions and verification boundary

An algorithmic claim needs more than existence of PSD blocks. The theorem
should quantify a finite-dimensional log-linear feature family, a supplied
bounded-width tree decomposition of the local-energy graph, a fixed bounded
parameter box (or explicit bit-cost/Lipschitz bounds), bounded local rates and
norms, and a certified exponential-oracle error budget. Under these
conditions, DP separation plus hypograph convex optimization can provide a
weak-optimization guarantee, with exponential dependence on treewidth and
local alphabet size.

The current numerical prototype does not establish that guarantee. A scipy
scalar optimizer or grid search is heuristic discovery, not a certified global
optimizer. Exact verification is separate: rational local ratios and DP
messages permit independent checking of every local PSD block and the summed
residual.

The positive control should use a known flat local-energy parameter so that
verification tests the bridge rather than optimizer luck. Evaluating a uniform
family at 39 bounded parameter points is a regression check, not a novelty
claim or a classifier for chemical Hamiltonians. Checkability alone cannot
imply chemical accuracy; that requires a separately proved approximation
condition for the chosen positive-factor family.
