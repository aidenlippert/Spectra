# Exact spin reorganization and the accepted H8 certificate

## 1. The claim that is actually accepted

For the supplied electronic Hamiltonian H on N=8 particles, the newly checked lower is

\[
L=-\frac{277669368076014353}{30000000000000000}\ {\rm Ha}.
\]

The upper is the exact expectation of the unchanged rational MPS, independently recomputed. The resulting width is 0.7674481354629528 mHa.

Acceptance uses the original symbolic certificate rather than assuming any numerical property of the new optimizer. The core identity is

\[
H_{\rm aux}-bI=\sum_\alpha B_\alpha^\dagger B_\alpha+(\widehat N-8)X+R.
\]

Here b=-9.255633788088 Ha after rational export. The exact spin-averaged coefficient remainder costs

\[
\eta=\frac{354426294353}{30000000000000000}\ {\rm Ha}.
\]

The new singlet lower is b-eta. The separate inherited M_S=1 bound is higher, so the new singlet controls the full fixed-N lower. The original-H spin defect, 59/250000000000 Ha, is charged once. Numerically proposed PSD matrices are exported as integer factor rows, making every accepted positive term an exact square.

No numerical dual optimum, excitation gap, discarded eigenvalue weight, or new checker rule is required for this conclusion. The certificate remains valid even if a proposal-space equivalence assertion were mistaken: the accepting checker independently expands the actual exported operators and charges their actual remainder.

## 2. The mixed cubic operator representation

Let J_+(B)=[S_+,B] be the adjoint spin action, with

\[
S_+=\sum_p a_{p\alpha}^\dagger a_{p\beta}.
\]

CAR gives

\[
[S_+,a_{p\beta}^\dagger]=a_{p\alpha}^\dagger,
\qquad [S_+,a_{p\alpha}]=-a_{p\beta},
\]

and zero for the other two spin choices. The action on a word follows the ordinary commutator product rule. Canonical reordering supplies only integer signs.

For each particle charge q=+1 or -1 and each of the two conserved spatial parity classes, the actual H8 mixed cubic dictionaries form four magnetic-component spaces. Ordered from m=-3/2 to m=+3/2, their dimensions are

\[
112,\quad368,\quad368,\quad112.
\]

All 3,840 mixed words across both charges were checked against literal CAR commutators. The integer raising/lowering maps also satisfy the spin commutator identity on all these spaces.

## 3. A sparse rational highest-weight basis

The raising map U from m=+1/2 to m=+3/2 is an integer 112-by-368 matrix. In these dictionaries,

\[
UU^T=3I_{112}.
\]

Its nonzero rows each have three signed entries and disjoint column supports. There are 32 columns not touched by U.

Every untouched column supplies a kernel vector. For a row supported at three columns with signs s0,s1,s2, use

\[
v_1=(s_1,-s_0,0)/2,\qquad
v_2=(s_0,s_1,-2s_2)/3.
\]

Both have zero row contraction and are orthogonal. These vectors, together with the 32 untouched columns, give 256 independent kernel columns K:

\[
UK=0,\qquad \dim\ker U=368-112=256.
\]

The implementation checks U(6K)=0 using integers. It checks that (6K)^T(6K) is diagonal with strictly positive entries. Thus no floating numerical rank threshold is involved in this kernel construction.

K spans the spin-1/2 highest-weight multiplicities. The m=+3/2 space supplies the 112 spin-3/2 highest-weight multiplicities. These are operator-coordinate spaces, not many-electron determinant spaces.

## 4. How an existing Gram construction is transported

Under exact SU(2) averaging, cross terms between inequivalent spin ranks vanish. The averages of equally normalized magnetic components in one irreducible multiplet agree. Thus positivity contributions can be represented using highest components, while retaining all Gram cross terms among multiplicities of the same spin.

For a quartet component, raise from m to +3/2. The squared normalization factors for m=-3/2,-1/2,+1/2,+3/2 are respectively

\[
36,\quad12,\quad3,\quad1.
\]

Applying the appropriate raising product on both sides of a covariance and dividing by these factors produces its highest-quartet covariance. For doublets, project the +1/2 component through K's exact left inverse; for the -1/2 component, first raise once, then project. The quartet contamination is removed by this projection.

Adding these transported covariances is positive. The old compact warm-start's complete projected coefficient vector is reproduced to 8.33e-17 maximum numerical discrepancy. Two random full-Gram tests additionally check the transport against the actual full coefficient maps. These floating map tests are construction diagnostics; the accepted energy comes from rational symbolic replay.

Unlike the preceding low-rank searches, this construction keeps the full 256-dimensional doublet and 112-dimensional quartet multiplicity spaces. It does not select their columns by MPS coverage or omit their internal cross terms.

## 5. Why the numerical system becomes sparse

Let M be the original spin-averaged physical word-pair coefficient map. With a sparse highest-weight basis K, the reduced Gram map is

\[
M(K\otimes K),
\]

symmetrized to respect the Frobenius inner product on real symmetric Gram matrices. The rational K has very few nonzeros. Consequently the construction need not materialize the dense projected coefficient tables produced by dense learned bases.

The complete new maps contain 1,370,610 nonzeros. If A is their assembled positive-cone map and F the free scalar/ideal map, the normal system is

\[
G=AA^T+FF^T.
\]

It has dimension 8,533 but only 221,779 stored nonzeros. Its sparse factorization has 954,952 nonzeros. A small regularization is used only in the numerical preconditioner; conjugate-gradient residuals and the true map are checked, and numerical error remains subject to exact final acceptance.

This is not automatically a lower asymptotic complexity theorem. The maps still derive from the full cubic word-pair preparation, and the number of retained multiplicity coordinates grows with system size. The immediate result is a successful sparse formulation of this particular H8 calculation.

## 6. The numerical state is a guide, not the physical constraint

The uploaded MPS was used to compute its actual degree-six moments. They initialize a proposal-side physical functional; no sextic moment was inferred solely from quartic moments. The lower optimization does not constrain the unknown ground state to equal the MPS.

The bound-maximization stage uses boundary-point updates in the sparse full-spin formulation. A subsequent lower-penalty refinement reduces reconstruction error. The accepted candidate was obtained from a checkpoint at 117.8004 s followed by 58.4562 s of refinement/export. Earlier experiments and original input discovery are separate costs.

The optimizer's final dual still does not provide a verified family optimum. That is unnecessary for a valid lower bound and forbids claiming that this energy is the optimal member of the cone.

## 7. Limitations and significance

The complete Gram-entry count is 353,296 including the separate nonsinglet proof. This is larger than the preceding failed compact construction's 115,032 entries. The new result therefore must not be described as a smaller version of that same search space. It trades a larger, complete spin-adapted space for sparse algebra and a successful bound.

The construction did not use the stronger full-cubic certificate's factors or eigendirections as a teacher. It did use the existing compact candidate as a warm start and inherits the MPS, nonsinglet proof, and original coefficient-map preparation.

The exact certificate closes H8's 1.6 mHa requirement with substantial margin. It does not establish the undisclosed Nooterra method, arbitrary-molecule scalability, experimental accuracy, or a theorem that the earlier compact spans could not have succeeded.

## 8. Final small-factor export

The initial accepted export contained 580 rows. Of these, 377 had integer coefficient L1 size below 1000 at denominator 10^9. The exact sum of their square operator-norm upper bounds is 198941/125000000000000000 Ha (approximately 1.59e-12 Ha). Their removal leaves 203 rows and 61,888 nonzero factor coefficients. The final checker re-expanded this smaller certificate and obtained the endpoints stated above. This reduces accepting representation size; it does not retroactively reduce the 335,168 singlet Gram entries used during numerical discovery.
