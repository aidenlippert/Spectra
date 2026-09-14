# Scope and soundness

Every reported lower endpoint is produced by the existing standard-library
CAR checker from an ordinary rational SOS certificate. The new numerical
search and moment code are untrusted proposers. All energies refer to the
supplied rational electronic Hamiltonian in the full fixed-N sector. The
MPS is an M_S=0 trial state; its moments never constrain unknown states.

## Complete coefficient reconstruction

Let w_i be one of the original charge-homogeneous operator dictionaries.
For a fixed real map V, write p_a=sum_i V_ia w_i. A PSD matrix S represents
sum_ab S_ab p_a^dagger p_b. The search includes every independent real
Hermitian CAR coefficient of degree 0, 2, 4 and 6 with the Hamiltonian's exact
charge/parity symmetries. The complementary conjugate coefficients follow
from real symmetry; the accepting checker expands and checks both. The
sector multiplier uses the same symmetry-filtered scalar, quadratic and
quartic Hermitian basis as the parent full-cubic problem.

The numerical equations are H=bI+SOS+(Nhat-N)X, with no free reconstruction
slack. After PSD clipping and rational factor rounding, the checker computes
the *actual* residual R. Since each CAR monomial is a contraction,
||R|| <= sum_w |R_w| = eta, and L=b-eta is a valid lower endpoint. A solver
status, a small numerical residual, or a low state action does not accept the
1.6 mHa requirement. It requires U_frozen-L <= 0.0016 exactly.

All degree-two blocks remain available as explicitly charged reconstruction
repair directions, including their cross terms. The teacher control picks
large eigenvectors of the accepted Gram matrices once per fixture. The
independent unpaired control picks small eigenvectors of each numerical MPS
moment matrix. The paired rule instead uses the sum of the two adjoint
moment matrices; it is a restricted subfamily of the unpaired cubic cone.
Its measured failure does not delimit the full cubic family.

## Paired cancellation and exact equality elimination

For odd p, p^dagger p + p p^dagger has no degree-six part when deg(p)<=3.
Mixed linear/cubic cross terms have degree at most four as well. When forming
the coefficient map, pair the indices of the adjoint block in transposed
order. The integer CAR maps then cancel before dense projection. The first
implementation omitted that transpose; its assertion refused the run. That
failed attempt is retained. Final integer exports explicitly permute the
same rounded factor, so pairing is exact after rounding too.

In this paired subfamily, reconstruction forces the degree-six part of
(Nhat-N)X_4 to vanish. `ideal_rank.py` generates that map from exact CAR
products. Gaussian elimination modulo the verified prime 65521 establishes
full column rank for the symmetry-restricted quartic multiplier basis.
A nonzero square minor modulo a prime is also nonzero over Q; therefore the
only quartic multiplier satisfying these zero equations is X_4=0.
This removes variables proved zero and identically zero equations, preserving
the paired model's exact, zero-remainder feasible set. The general certificate
format allows a nonzero remainder and quartic X; the equivalence does not
extend to every such residual-corrected certificate. It also makes no such
claim for unpaired blocks. A four-mode counterexample tests refusal when the rank
argument fails. This is an exact algebraic construction, not a theorem that
the compressed proof family cannot meet an energy target.

## What the cut construction implements

Within each odd dictionary, partition words by their prefix to the left of
two orbital cuts (at one third and two thirds of that dictionary's orbital
extent). Starting from the finer cut, retain at most two low-action
combinations per prefix node, then merge at the coarser cut. Finally retain
the requested root rank. Child maps have disjoint supports and are checked
as isometries. Their composition gives a fixed map on the original CAR
operators. The root rank and per-node rank are independent of MPS bond size.

For every physical state rho, C_rho[i,j]=Tr(rho w_i^dagger w_j) is PSD.
Consequently V^dagger C_rho V is PSD for every fixed composed map, including
one chosen using an entirely different trial state. Only the root congruences
are retained as PSD constraints; intermediate nodes construct the maps and
are not additional PSD constraints. The root constraints share the same
global CAR coefficient rows and sector relations.
No product-state, nearest-neighbor, fixed-density, or MPS restriction is
imposed on rho. The full long-range Hamiltonian remains in the energy row.

This implementation is a finite two-level operator-map restriction. It still
builds the full dictionary moment matrices and original sparse CAR maps;
its dual export is expanded into ordinary factors. It does **not** implement
an efficient tensor contraction of a molecular RDM consistency hierarchy, or
demonstrate that such a hierarchy avoids full coefficient-map construction.
The cut-map experiments must be assessed with this limitation and their
preparation/replay costs included.

The relevant precedent is Kull et al., [Lower bounds on ground-state energies
of local Hamiltonians through the renormalization group](https://arxiv.org/abs/2212.03014),
which develops compressed consistency constraints and discusses an RDM
extension. [Coarse-grained Bootstrap of Quantum Many-body Systems](https://arxiv.org/abs/2412.07837)
also uses tensor maps for spin-chain constraints. Neither paper supplies the
finite molecular implementation claimed or tested here; tensor coarse
graining itself is not being claimed as new.

## Rigorous upper verification

The new checker reuses every structural/charge-flow gate of the integer MPS
checker, then establishes norm positivity using its own rigorous contraction.
Each exact rational tensor and Hamiltonian coefficient is enclosed using an
exact Fraction comparison with the converted float. Each binary64 addition
and multiplication is rounded outward to adjacent floating-point values.
Sign changes and exact zeros are handled explicitly. Shared-prefix transfer
contractions are evaluated with fast-math and reassociation disabled, and
nonfinite enclosures are refused. IEEE-754 binary64 with gradual underflow
is the stated arithmetic requirement. This is an interval proof, not an
ordinary floating-point energy estimate.

Acceptance requires a strictly positive lower norm enclosure and a
nonnegative lower enclosure for u*norm-<phi|H|phi>. An ambiguous sign refuses
or invokes the existing exact integer checker. In particular an endpoint
exactly equal to the rational Rayleigh quotient usually requires fallback.

For a matched speed comparison, both checkers receive the same state and
u=ceil(2^24 U)/2^24. The added allowance is explicit and below 2^-24 Ha.
The main compressed-lower comparisons continue to use the original exact U.
The integer oracle separately verifies enclosure containment and the same
comparison endpoint. Tests include invalid endpoints, a zero state, malformed
charge flow, altered bindings, duplicate entries, denominator errors, and
primitive intervals compared with exact rational arithmetic.

The arithmetic APIs and compiler assumptions are documented in the official
[Python math documentation](https://docs.python.org/3/library/math.html#math.nextafter)
and [Numba floating-point semantics](https://numba.readthedocs.io/en/stable/reference/fpsemantics.html).
