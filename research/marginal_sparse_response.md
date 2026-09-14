# Sparse molecular response with an exact residual certificate

The H4 square's sparse-construction precision gap is closed on the stored rational fixture. Keeping the same P32 reference, occupation-tree complement proof, and integer upper witness, **12 response directions** improve the certified interval from **0.16418909584125116 to 1.002511669462418e-10 Hartree**.

Construction and replay use sparse Hamiltonian actions and small Gram matrices. They do not form or invert the complete Q matrix or diagonalize the complete sector. However, final replay evaluates **all 70 determinant source actions**. This combines sparse construction with directional response; it does not establish a compact configuration frontier or favorable molecular scaling.

## Residual identity

Let the already checked occupation tree certify H_QQ>=gamma I. For a candidate lower endpoint b<gamma, define

\[
A=H_{QQ}-bI,\quad\delta=\gamma-b>0,\quad W=H_{QP}.
\]

For any Q-supported response approximation X, write R=W-AX. The exact identity

\[
W^TA^{-1}W
=W^TX+X^TW-X^TAX+R^TA^{-1}R
\]

and A^(-1)<=delta^(-1)I imply the upper response

\[
\Sigma(X)=W^TX+X^TW-X^TAX+\delta^{-1}R^TR.
\]

The remaining response is enclosed by a rigorous residual term. X need not be an exact inverse action or lie in an invariant subspace.

## Optimize the bound in a small response space

Choose an exact sparse basis K and put X=KZ. Let B=A-delta I=H_QQ-gamma I>=0. Since A and B commute, AB>=0. Define

\[
D=K^TABK,\qquad L=K^TBW.
\]

Expanding the residual gives

\[
\Sigma(KZ)=\delta^{-1}
\left(W^TW-L^TZ-Z^TL+Z^TDZ\right).
\]

When D is positive definite, completing the square gives the Loewner-minimal choice Z=D^(-1)L and

\[
\boxed{\Sigma_K=\frac{W^TW-L^TD^{-1}L}{\gamma-b}}.
\]

Thus the exact sufficient lower certificate is

\[
H_{PP}-bI-\Sigma_K\succ0.
\]

Both D and this retained matrix are checked using exact rational positivity. The implementation rejects a singular D instead of applying an unchecked pseudoinverse. Full rank of K is checked separately through K^T K. The complete complement proof remains mandatory, including uncoupled Q states.

For fixed gamma and b, this optimized bound is never worse than the original scalar W^T W/(gamma-b). Enlarging the response span can only decrease the optimized upper response. If the span contains A^(-1)W, the residual vanishes and the exact response is recovered. These are algebraic statements about the bound, not guarantees that a small useful span will be discovered efficiently.

## Only action Grams are required

The exact preparation computes

\[
G=K^TK,\quad T=K^TH_{QQ}K,\quad
V=(H_{QQ}K)^T(H_{QQ}K),
\]

\[
C=K^TW,\qquad J=(H_{QQ}K)^TW.
\]

Then

\[
D=V-(b+\gamma)T+b\gamma G,\qquad L=J-\gamma C.
\]

No inverse or spectral decomposition of H_QQ appears. Sparse integer response directions are serialized, and their actions and all Grams are reconstructed from the Hamiltonian during replay. Directions containing P support, duplicate states, zero vectors, invalid occupations, or dependent columns are rejected.

## Discovery and the controlled result

The proposer finds a violating retained direction in the current approximate Schur matrix. It applies W to that direction and approximates its inverse response by conjugate gradients using sparse Q actions. The candidate is numerically orthogonalized, normalized, and rounded into a sparse integer vector. Exact replay then checks the rounded vector and the entire certificate. Neither the numerical iteration's residual nor its convergence flag is accepted as a proof.

The response budget is 16 directions; this run succeeds at 12. Its response support sizes are 4,4,12,9,12,6,6,14,14,4,4,6. The retained matrix stays dimension32 and the only exact response solve is dimension12. The 38-dimensional Q inverse is never formed.

| Response dimension | Certified width (Ha) |
|---|---:|
| 0, original scalar bound | 0.16418909584125116 |
| 1–9 | 0.16418909584125116 |
| 10 | 0.16418909512125116 |
| 11 | 0.035857644396251165 |
| 12 | 1.002511669462418e-10 |

The unchanged early widths are recorded, not discarded. Several unresolved directions must be addressed before the lower endpoint improves substantially. This is not evidence that single-direction selection always converges quickly.

The certificate is 55,907 bytes. Initial proposal plus exact checks took about 2.13 seconds; the recorded standard-library replay took about 0.10 seconds. These are finite-run timings. The 12-direction proof reaches the same approximately1e-10-Hartree margin used by the conservative lower proposer; it is not a claim of physically meaningful accuracy beyond the exported rational Hamiltonian.

## State coverage and cost accounting

The proposal's own cache evaluates and references 68 states. The final independent replay additionally performs the occupation-tree leaf checks and evaluates all70 states. Construction totals include those replay checks: **70 source states and 70 referenced determinants**, not just the smaller proposal cache. The original sparse-construction square already referenced all70 destinations, although it evaluated only65 source actions.

The method therefore removes the full matrix/inverse requirement while retaining a full source frontier on this test. A sparse data structure is not itself a complexity result. Dense response Grams, retained dimension, rational coefficient size, and the action frontier all remain relevant costs.

## Verification and unresolved transfer

An independent `python -S` replay accepts the export at `results/marginal_sparse_molecular/square_residual_response/certificate.json`. Tests establish, on a noncommuting rational example, both the upper-response inequality and monotonic improvement with nested response spaces; the full response space reproduces the exact inverse. Another test directly expands the residual identity and compares it with the optimized Gram expression.

Rejection tests cover nonpositive gaps, independent bases with singular D, dependent response vectors, P-contaminated support, missing response recipes, false lower endpoints, and incomplete occupation coverage. The production proposer is rerun with the full-matrix builder forbidden, retained eigensolves limited to32, and response solves limited to16. Replay is also checked with its numerical response proposer disabled.

The next unresolved issue is simultaneous control of the reference, response, and state frontier as the interacting system grows. This fixture preserves the preceding upper witness; a perturbation that puts appreciable ground-state weight outside P will also need a strengthened sparse upper witness. The subsequent [H6 transfer](marginal_h6_transfer.md) proves that such an upper adaptation is necessary for its selected P32, supplies sparse physical Krylov uppers, and isolates a complement row-bound failure on924 states. The subsequent [H6 complement proof](marginal_h6_complement.md) proves that no positive diagonal weighting can repair that threshold, replaces it with a complete sign-preserving factor certificate, and reaches width1.22199e-10 Ha with32 response directions. An exact fixed-target rank obstruction requires at least31 directions in that response family. Full-sector coverage remains; uniform error/cost bounds, the asymmetric quartic optimum, and a universal physical-marginal boundary remain unproved.
