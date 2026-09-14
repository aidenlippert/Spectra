# A gapped vacuum does not constrain a restricted pseudomoment dual

For a path of M fermionic modes, let

H = sum_i n_i - t sum_i (a_i^dagger a_(i+1) + h.c.),  0 < t < 1/2.

The one-particle matrix is I minus t times the path adjacency matrix, whose
spectrum lies in [-2,2]. Every one-particle energy is at least 1-2t. On the
full Fock space the unique ground state is the vacuum, with a uniform gap at
least 1-2t and zero connected ground-state correlations on disjoint supports.
These are analytic statements; no many-body diagonalization is used.

Consider the lower-bound relaxation H = bI + sum_w lambda_w w^dagger w + R,
with lambda_w >= 0, coefficient-L1 residual penalty, and only SINGLE MONOMIAL
squares retained. No particle-number ideal is imposed in this example.
Its dual requires y(I)=1, |y(monomial)|<=1, y(w^dagger w)>=0 for retained w.
It does not require positivity on arbitrary linear combinations of words.

Set y(I)=1, all nonconstant diagonal occupation moments to zero, and
 y(a_i^dagger a_j)=1 for nearest neighbors and for endpoints (0,M-1),
with their Hermitian partners. All other normal-ordered moments are zero.
Every monomial square is diagonal in occupations, so its y value is precisely
its nonnegative vacuum value. The functional is feasible for the restricted
cone, but is intentionally nonphysical.

It has y(H)=-2t(M-1). This is the restricted optimum: the primal b=0,
SOS=sum_i n_i has residual equal to the hopping terms, with coefficient-L1
norm 2t(M-1), and therefore attains the same objective.

The omitted endpoint block for words (a_0,a_(M-1)) is [[0,1],[1,0]]. Its
normalized Rayleigh quotient for (1,-1) is -1, independent of the distance
M-1. Thus a uniform gap and even exactly vanishing physical correlations do
not imply distance-decaying pricing violations in a restricted dual.
This does NOT prove an objective gain from adding that endpoint block, nor
refute a locality theorem for a stronger, explicitly constrained dual cone.
In particular it does not establish failure for the current number-ideal
molecular relaxation; it refutes using the physical gap alone as its proof.

The root-corrected stdlib script uses the actual CAR word_product, checks all
normal-ordered monomial squares through degree two (37,137,301,529 at
M=4,8,12,16), computes the endpoint block, and matches exact primal/dual
objectives. Receipt: results/certificate_scaling/locality_dual_exact.json.
Earlier locality_dual_counterexample.json and agent-only concatenation checks
are superseded and must not be used as evidence.

An orbital graph is also different from an operator-word interaction graph.
Sparse quantum SOS requires a proof on the latter, including overlap and
quotient constraints. Neither a small one-particle treewidth nor a small
physical correlation length supplies that proof automatically.
