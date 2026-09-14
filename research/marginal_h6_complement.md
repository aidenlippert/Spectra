# H6: comparison obstruction, exact complement gap, and global interval

The rational H6 electronic Hamiltonian now has a certified ground-energy interval of width **1.221992366627265e-10 Hartree**. The lower endpoint combines a sign-preserving complement factor proof with 32 optimized response directions. The upper is the previous independently evaluated sparse physical Krylov witness. This resolves the finite H6 precision blocker; that certificate explicitly evaluates all 924 fixed-number determinants. The subsequent [sharper threshold and spin reduction](marginal_h6_spin_reduction.md) reduces the successful response to28 directions and final replay to400 distinct configurations, with an original-H interval width1.94199e-10 Ha. Discovery still inherits the full-sector proof; efficient scaling remains unproved.

## The absolute-value approximation is provably insufficient

Keep the previously selected coordinate reference P of dimension 32. Let Q be its 892-dimensional orthogonal complement. For the real symmetric matrix H_QQ, define its comparison matrix by

\[
C_{ii}=(H_{QQ})_{ii},\qquad C_{ij}=-|(H_{QQ})_{ij}|\quad(i\ne j).
\]

For every strictly positive diagonal weighting d, the weighted Gershgorin lower bound is

\[
g(d)=\min_i (Cd)_i/d_i\le\lambda_{\min}(C).
\]

This follows by applying Gershgorin to the diagonally similar matrix diag(d)^{-1} C diag(d). A nonzero vector x therefore gives an obstruction whenever

\[
g(d)\le\lambda_{\min}(C)\le x^TCx/(x^Tx)<\gamma.
\]

The exported integer witness has only 64 nonzero coordinates. Exact CAR action and rational Rayleigh evaluation place it **0.4804954363965175 Ha below** the requested threshold

\[
\gamma=-6.323058626212\ \mathrm{Ha}.
\]

Consequently, **no positive diagonal weighting of determinant-basis Gershgorin can certify this threshold for this P**. The obstruction permits disconnected Q and zero witness amplitudes; it requires no Perron-vector reconstruction. It concerns the comparison family, not the actual physical spectral gap.

Separate numerical controls on explicit 892-state sparse matrices estimate the physical Q minimum at -6.2647991010 and the comparison minimum at -7.9852472113 Ha. These eigenvalues are diagnostic, not certified endpoints. They identify sign cancellation as a major loss in the comparison bound.

## Certifying the physical complement while retaining signs

The exact Q action graph splits into blocks of dimensions

\[
117,200,108,168,117,108,18,18,18,18,1,1.
\]

Replay validates every state, excludes P, checks uniqueness and the exact cardinality binomial(12,6)-32, and rejects any nonzero interblock Q coupling. This proves full complement coverage, including blocks with no coupling to P. Every Q coordinate is explicitly listed.

For each block, put A=H_block-gamma I. The integer factor certificate represents a unit lower triangular L and positive diagonal D at a common scale S. Its entries are L_ij=l_ij/S below the diagonal, L_ii=1, and D_ii=d_i/S. A itself must be represented exactly at scale S; replay never rounds the operator.

With implicit l_ii=S and Abar=SA, replay reconstructs the residual using integers:

\[
\overline E_{ij}=\overline A_{ij}S^2-\sum_k l_{ik}d_kl_{jk},
\qquad e=\max_i\sum_j|\overline E_{ij}|/S^3\ge\|A-LDL^T\|_2.
\]

Upper bounds on inverse triangular row and column sums follow from

\[
r_i=S+\left\lceil\sum_{j<i}|l_{ij}|r_j/S\right\rceil,
\quad
c_i=S+\left\lceil\sum_{j>i}|l_{ji}|c_j/S\right\rceil.
\]

Thus ||L^{-1}||_2^2 <= r_max c_max/S^2. The accepted positivity margin is

\[
\mu=\frac{d_{\min}S}{r_{\max}c_{\max}}-e>0.
\]

All 12 blocks pass at S=10^12. The smallest verified margin is approximately 2.56894e-5 Ha. This certifies H_QQ strictly above gamma without treating rounded factorization as an exact identity. Proposal rounding uses signed nearest-ties-even integer division; acceptance recomputes the entire residual independently.

The standalone gap replay evaluates 892 source actions and references all 924 determinants. Factor construction plus the separate comparison diagnostics initially took about 3.86 seconds. These are finite control costs, not a scaling result.

## The resulting global energy proof

Use A=H_QQ-bI, delta=gamma-b>0, W=H_QP, and an integer response basis K contained in Q. The existing optimized residual bound gives

\[
\Sigma_K=\frac{W^TW-L_K^TD_K^{-1}L_K}{\delta}\succeq W^TA^{-1}W,
\]

where

\[
D_K=K^TA(A-\delta I)K,\qquad L_K=K^T(A-\delta I)W.
\]

Exact sparse actions reconstruct the Gram matrices. Replay checks independence of K, strict positivity of D_K, and strict positivity of H_PP-bI-Sigma_K. It uses the same decoded Hamiltonian, P, and gamma for the complete complement factor proof. A separate exact Rayleigh evaluation supplies u.

The final electronic energy endpoints are

\[
b=-6.333058626334,
\qquad
u=-6.3330586262118007633\ldots\ \mathrm{Ha}.
\]

The exact fractional upper and width are stored in the receipt. Unlike agreement between two upper witnesses, these endpoints bound the actual ground energy of the exported rational Hamiltonian over the full fixed-number sector.

| Response dimension | Certified interval width (Ha) |
|---:|---:|
| 16 | 0.10682918096919923 |
| 24 | 0.06452648454419924 |
| 30 | 0.025953919469199237 |
| 31 | 0.019535073952199236 |
| 32 | 1.221992366627265e-10 |

The response vectors have 160 or 168 nonzero Q amplitudes. Complete replay evaluates and references all 924 determinants. The certificate occupies 1,938,902 bytes. The first 16-direction pass took about 74.00 seconds; continuation through dimension 32 took about 209.59 seconds. These phase times include repeated exact acceptance checks. The numerical response proposer uses sparse conjugate-gradient actions; exact replay invokes neither CG nor an eigensolver.

## A precise response-size obstruction

At the fixed target b_star=u-10^-10 and the certified gamma above, define the scalar Schur bound

\[
S_0=H_{PP}-b_\star I-W^TW/(\gamma-b_\star).
\]

Deleting retained-coordinate index 8 leaves a 31-by-31 principal matrix whose negative is exactly positive definite. Therefore S_0 has at least 31 negative directions. The improvement to S_0 supplied by r response directions is PSD and has rank at most r. The intersection of this negative subspace with the correction kernel remains nonzero when r<31. Hence **fewer than 31 response directions cannot certify this target with this scalar-plus-rank-update formula**. This is a family-specific exact lower bound, not a universal lower bound on chemical certificate size.

The complementary ideal upper bound is also clear. If K contains range(A^{-1}W), the residual response is exact. This space has dimension rank(W)<=32. At a fixed target with positive exact Schur matrix, targeting a nonpositive direction of the current bound must add an independent exact inverse-response vector: if it were already in K, the bound would already agree with the positive exact Schur form on that direction. Ideal exact targeting therefore terminates within rank(W) updates. Approximate solves and integer rounding require the implemented replay gates and do not inherit exact finite termination automatically.

## What remains

The finite H6 lower-energy and precision blockers are resolved. The construction still pays for the entire Q sector and factors blocks as large as 200. Both eliminating that complete enumeration and controlling response dimension as system size grows remain open.

The numerical sensitivity probe identified -6.265 Ha as a sharper threshold. The subsequent compression experiment now certifies that threshold and reaches width1.40199e-10 Ha with28 response directions. An exact nearby-SU(2) reduction then lowers final replay to400 distinct configurations, with width1.94199e-10 Ha after perturbation transfer. The numerical negative count of21 is not yet an exact rank certificate, and no21-direction interval is claimed. The subsequent [direct spin constructor](marginal_h6_direct_spin.md) now removes the inherited full-sector discovery cost and certifies width1.93199e-10 Ha from the Hamiltonian alone. A perturbation merging the Q blocks is the next concrete transfer test; its reference gap survives an exact norm shift.

The result concerns a rational finite-basis electronic Hamiltonian. It does not certify integral computation, basis completeness, finite-temperature properties, synthesis pathways, or a universal boundary representation for physical marginals. Explicit sector factors have not removed the many-body complexity.

## Replay

```sh
python -S -m experiments.marginal_h6_complement --verify-comparison results/marginal_h6/complement_proofs/comparison_certificate.json
python -S -m experiments.marginal_h6_complement --verify-factor results/marginal_h6/complement_proofs/factor_certificate.json
python -S -m experiments.marginal_sparse_response --verify-rank-obstruction results/marginal_h6/response_rank_obstruction/certificate.json
python -S -m experiments.marginal_sparse_response --verify results/marginal_h6/complement_proofs_response_continued/certificate.json
```

All four replays passed without site packages. The intermediate 16-direction interval also passed independent exact replay. Focused tests independently reconstruct factor residuals and inverse triangular norms, exercise signed rounding and false-factor rejection, require complete Q coverage and zero interblock coupling, check physical-Hamiltonian binding, validate the rank obstruction, reproduce a small factor-response construction and continuation, and replay the final H6 interval.
