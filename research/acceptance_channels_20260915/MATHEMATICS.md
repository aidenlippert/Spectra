# Acceptance objective and collective paired channels

The certified object is the same exported rational H12 Hamiltonian. All statements
below concern that finite model and the explicitly specified spin/number domain.

## What the H12 receipt establishes

The previous full original-model width was 12.590172568672 mHa. Its singlet
coefficient remainder allowance was 10.648044614770 mHa (84.5743%). Subtracting
only that allowance *hypothetically* leaves 1.942127953902 mHa. This subtraction
is not a new certificate. The residual is part of the exact identity.

The old numerical solver already chooses its export by b minus the accepting
coefficient-L1 allowance. Its iterations optimize an equality-constrained SDP,
and its ideal recovery uses a least-squares quotient. The new LP experiment
changes the ideal fit itself, rather than relabeling the old export score.

## An explicit accepting-objective fit

Write the row-scaled reconstruction as A(Q)+F x=h. The first coordinate of x is b.
Let T be the stored spin projector, S its independent selected rows, and D the
row-scaling matrix. The unique invariant residual lift is

    K = T[:, S] (T[S, S])^{-1} D^{-1}.

The numerical accepting objective is

    b - sum_i w_i |[K(h-A(Q)-F x)]_i| - spin_allowance.

The weights restore the omitted Hermitian partner of each representative row.
An exact export and full original-model replay are still required.

The fixed-Q LP optimizes b and the independent ideal coordinates. The block-scale
LP also replaces Q_j by a_j Q_j with a_j >= 0. Both are convex restrictions of
the existing full Gram problem. Epigraph variables encode each absolute value.
A numerical optimum of either LP is not an exact full-family obstruction.

For H12 the residual lift consists of 37,566 disconnected spin components, largest
dimension five, with six distinct inverse matrices. This is a direct reuse of
the existing projector's structure; it does not prove cheap global certification.

The full SCS experiment removes the restricted-subcone assumption. Using

    ||K r||_(w,1) = max_(|d_i|<=w_i) d^T K r,

its conic dual has variables y and d and the constraints

    minimize h^T y
    subject to F^T y=e_b, A^*(y)>=0, y=K^T d, |d_i|<=w_i.

These are the complete original Gram cones plus explicit residual-box variables.
The numerical program has 414,768 scalar variables, 1,389,349 conic rows and
3,586,991 nonzeros on H12. These extra sizes and their costs are charged; the
acceptance-aware reformulation is not computationally free. The conic dual
multipliers recover all Gram matrices and ideal coordinates. PSD projection,
recomputed residual scoring, rational export and the original exact check are
still required. A time-limited SCS output is not an exact optimality certificate.

## The pasted H8 channel in the current orbital coordinates

The supplied B satisfies B=A a_(3 beta), and

    B†B + BB† = A A† + n_(3 beta) [A†, A].

Its exact expansion has 31 terms and degree at most four. The current H8 and H12
highest-weight dictionaries contain B. They also contain the spin-rotated B†,
called C. Exact spin averaging gives

    T(B†B + BB†) = T(B†B + C†C).

The retained bases of these two blocks are exactly identity matrices. Explicit
rational coefficient vectors v and w therefore supply the PSD Gram additions
v v^T and w w^T. The executable receipt `channel_membership.json` gives all
coordinates and source hashes. Consequently this isolated ray cannot enlarge
either *current* cone. The older external paired-window cone is different.
Its source dual was not supplied, so its claimed negative evaluation is not
independently replayed here.

This membership statement interprets the supplied polynomial in the current
dictionary's orbital coordinates. The external bundle did not identify an
accessible source Hamiltonian or orbital transformation here. If its original
operator uses a different orbital basis, that physical operator must first be
transported; membership of the transported operator has not been established.
The H12 channel searches use their own frame-bound moments and do not depend on
this external basis identification.

## A larger but compact set of collective coordinates

Use alpha modes a_p and beta modes b_q. For real coefficients c_(p,q,r), q<r,
define

    B = sum_(p,q<r) c_(p,q,r) a_p† b_q b_r.

There are 792 such monomials for twelve spatial orbitals. A paired direction
B†B+BB† is positive on every physical state. Its sextic terms cancel exactly.
No sixth-order moments are needed to test this paired condition.

Set C_(p,(q,r))=c_(p,q,r), and extend it antisymmetrically as
A_(p,q,r)=-A_(p,r,q), with zero diagonal. Then the exact quartic expansion is

    sum_(p,s) (C C^T)_(p,s) a_p† a_s
  - sum_(q<r,t<u) (C^T C)_((q,r),(t,u)) b_t† b_u† b_q b_r
  - sum_(p,s,v,w,k) A_(p,w,k) A_(s,v,k) a_p† b_v† b_w a_s.

For cross terms between two collective directions, symmetrize the corresponding
bilinear contractions. These are Gram products of matrices with sizes 12x66,
66x12 and 144x12, rather than a materialized global cubic coefficient map.
Integer numerators are contracted with an explicit overflow check; floating
division is proposal-only. Exact polynomial tests compare the contractions with
literal CAR multiplication, including all cross terms.

The selected directions are columns of V. A complete small Gram Q >= 0 yields
factors from sqrt(Q) V^T, with a matched adjoint factor block. Thus the directions
can interfere through off-diagonal entries of Q. Four directions add 16 Gram
entries, but each direction may contain 792 coefficients. These are separate
representation-size measures.

The reference accepting path remains the original standard-library checker. It
re-expands the exported integer factors; the contracted numerical map cannot
substitute for that calculation. The separate exact accepting-cost prototype
below uses integer contractions and requires a full reference comparison.

## Normal-system update without a dense outer product

If M contains the new projected Gram columns, the numerical normal operator is

    G_new = G_old + M M^T.

The implementation applies this expression as an operator. Its preconditioner
uses the Woodbury formula on the already required regularized G_old factorization:

    (G_old + M M^T)^-1
      = G_old^-1 - W (I + M^T W)^-1 M^T G_old^-1,
    W = G_old^-1 M.

The formula is tested against an explicitly formed SPD matrix. The existing
normal-map consistency check and iterative solve tolerances remain in force.
This is a computational representation of the same new normal equations, not
an approximation used by the exact verifier.

## Scope of the selector

The source dual is approximate. Its positivity violations propose directions;
they do not prove the existing family's optimum. A failed truncated-vector
search likewise proves no sparsity obstruction. The 792x792 moment matrix is
charged discovery work. Reusing a handful of its directions is not evidence of
a matrix-free selector or general scaling.

The relevant success criterion is a tighter independently accepted full interval,
with selection, preparation, failed searches and replay costs reported.

The second spin channel uses general modes c_i and words c_p† c_q c_r with q<r
and magnetic-number change one. The same three contractions remain valid with
the indices running over all 24 spin orbitals. Coinciding creator/annihilator
indices are allowed; repeated creators or repeated annihilators in a resulting
normal word vanish by CAR. Exact tests include those overlaps.

Its diagnostic matrix has 2,520 rows. Every entry comes from

    {c_p† c_q c_r, c_u† c_t† c_s}
      = delta_(p,s) c_u† c_t† c_q c_r
      + (delta_(r,u) delta_(q,t)-delta_(q,u) delta_(r,t)) c_p† c_s
      - delta_(r,u) c_p† c_t† c_q c_s
      + delta_(q,u) c_p† c_t† c_r c_s
      + delta_(r,t) c_p† c_u† c_q c_s
      - delta_(q,t) c_p† c_u† c_r c_s.

The vectorized moment matrix is tested against exact singlet traces of every
entry on a two-spatial-orbital case. A required missing moment causes rejection;
it is never filled by a guess. Four such directions accompany the prior
eight-direction block as a separate paired Gram block, preserving its complete
cross terms. The two blocks add 64+16=80 Gram entries. Their physical coefficient
lists, large selector matrices and unchanged base maps remain additional costs.

## Exact accepting-cost prototype

For an exported row of integer factor coefficients, the same three contractions
can be evaluated with unbounded Python integers before dividing by the squared
factor denominator. This yields exactly B†B+BB† through degree four. The
prototype first applies the original validation to every dictionary and row.
It requires the following dictionary to be the literal adjoint in the same
order and to have exactly identical factor rows. Unmatched or independently
weighted adjoints take the original expansion path.

Exact tests compare the whole polynomial, including its zero sextic remainder,
and the factor statistics against the original CAR implementation. They cover
multiple dense factor rows, duplicate and reversed words, arbitrary-size
integers, unequal adjoints and original refusal paths. A full molecular receipt
comparison is an additional requirement for the measured prototype. This
changes representation of an exact sum, not the allowed residual or any energy
acceptance condition.

## What the new exact receipt says about residual-only repair

For the replayed four-direction export, the full width is 2.588188188369 mHa
and the singlet remainder allowance is eta = 0.365389793467126 mHa. This is
14.12% of the width, rather than the old receipt's 84.57%.

Keep this exported b, its positive squares, the upper and all other allowances
fixed. Suppose only the residual lower -eta is replaced by a better scalar r.
The already checked norm bound implies -eta I <= R <= eta I, so every valid
scalar lower satisfies r <= eta. The possible improvement is at most 2 eta.
Even this optimistic substitution leaves at least 1.857408601435 mHa.

Thus a better scalar bound on this fixed remainder alone cannot close the
1.6 mHa target. Reoptimization of b and the squares, or a bound on their joint
spectrum with the remainder, is still allowed. This is not a family obstruction
and does not distinguish inadequate convergence from inadequate representation.

## Exact physical null channels behind part of the trace diagnostic

Let S_plus = sum_q a_(q alpha)† a_(q beta). Exact CAR gives

    [S_plus, a_(p beta)] = 0,
    [S_plus, a_(p alpha)†] = 0.

Since S_plus annihilates every singlet, the cubic operators

    Z_p^- = S_plus a_(p beta),
    Z_p^+ = S_plus a_(p alpha)†

also annihilate every singlet. Consequently their coordinate vectors are null
vectors of every physical singlet Gram matrix in a containing dictionary.
For the current H12 blocks 22 and 24, all twelve vectors in each block are
contained in the retained identity basis. Their exact coordinate Gram matrix
is 11 I_12, proving independence. The receipt saves all 24 rational vectors.
An independent ladder-action test checks the full three-dimensional singlet
basis for two electrons on two spatial orbitals; it also checks that these
operators do not vanish on every spin sector.

This proves physical null directions, not a repaired moment functional. It
does not establish which null equations are forced by the current finite ideal
map. Nor does it repair the other trace-kernel blocks. A nonzero action on the
physical trace kernel alone is not, without those extra conditions, a general
test of infeasibility of the original SDP dual. The proposed trace-mixing repair
still lacks exact affine/nullspace handling and accepted PSD checks.
