# A solvable six-mode stress test for sparse certificates

Consider fixed particle number `N=3` and

`H(t)=C(N_A,2)+C(N_B,2)-t sum_{i=0}^2(a_i†a_{i+3}+a_{i+3}†a_i)`,

where `A={0,1,2}`, `B={3,4,5}`. Direct exact Fock diagonalization gives

`E0(t)=2-t-sqrt(1+2t+4t^2)` for `t >= 0`.

This follows analytically in the fully permutation-symmetric sector. Map each
matched pair `(i,i+3)` to a flavor spin and let `Sx` be the collective flip
operator. The interaction is diagonal in `m=(N_A-N_B)/2`, with values 3,1,1,3
for `N_A=0,1,2,3`; hopping is `-2t Sx`. The lowest eigenvalue of the relevant
spin-3/2 block is the expression above. The remaining permutation sectors have
larger eigenvalues for `t>=0` (this comparison is a finite 20-dimensional
calculation and should be retained as a regression check).

The exact finite-sector certificate is therefore a block certificate: decompose by
the simultaneous `S_3` representation on the three matched flavors, then
factor each finite block `H_lambda-E0 I` as `R_lambda†R_lambda`. This is an
exact certificate on the finite sector, but there is no general low-degree
guarantee: pulling the factors back can require many-body projectors or
polynomials in conserved-sector projectors. In particular, at `t=0` the
degenerate occupancy boundary is a warning against claiming that DQG alone is
exact. Discovering a factorization whose summands reduce to quadratic words is
an additional, testable property—not a consequence of block diagonalization.
The square-root coefficient can be handled as an algebraic number
(`x^2=1+2t+4t^2`) or with interval arithmetic.

This is a deliberately restricted integrable family, not evidence of a generic
short certificate. Its value is diagnostic: it tests whether a compiler can
discover symmetry-adapted local occupancy blocks and produce a verifiable dual
factorization whose lower bound meets the explicit primal eigenvalue. For a
generic hopping matrix the simultaneous permutation reduction disappears; the
same block/SOS recipe is then a search heuristic rather than a closed form.
