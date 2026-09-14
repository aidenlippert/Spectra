# V5 deflation proposal: what can actually be claimed

Consider a three-dimensional Gaussian source with covariance

\[
C=\lambda_1u_1u_1^T+\lambda_2u_2u_2^T+\lambda_3u_3u_3^T+\sigma^2I,
\qquad \lambda_1>\lambda_2>\lambda_3>0,
\]

where the orthonormal frame (u_i) is selected from a declared finite rational set. An action chooses a unit vector (x), and the scalar observation has variance (v(x)=x^TCx), with independent Gaussian shot noise or a bounded sub-Gaussian estimator. No vector readout is available: every direction costs samples.

## Finite acquisition construction

Let the candidate frame set have minimum separation (ho), and let the eigengaps be (Delta_{12}=lambda_1-lambda_2) and (Delta_{23}=lambda_2-lambda_3). A first-stage action design estimates (u_1) by comparing candidate directions and requires enough samples that the variance-estimate error is below a gap-dependent margin (g_1). If (hat u_1) has angle (	heta) from (u_1), a nominally null probe (x\perp\hat u_1) retains leakage

\[
(x^Tu_1)^2\lambda_1\leq\lambda_1\sin^2\theta.
\]

Therefore deflation is useful only if (lambda_1sin^2\theta) is smaller than the second-mode action gap (g_2). Under this condition, the second-stage candidate set can be restricted to frames orthogonal to the certified (u_1) neighborhood, reducing the number of pairwise comparisons from (|\mathcal F|) to the surviving subset (|\mathcal F_1|). After (u_2) is certified, the same argument applies to (u_3), with residual leakage bounded by (lambda_1sin^2\theta_1+lambda_2sin^2\theta_2).

For sub-Gaussian variance estimates with proxy (R), a sufficient sample count for a finite candidate comparison is of order

\[
n_j\geq \frac{2R^2}{g_j^2}\log\frac{2|\mathcal F_j|}{\delta_j},
\]

where (g_j) is the smallest declared variance gap after accounting for prior leakage. This is a finite upper bound, not a dimension lower bound. Exact Gaussian tail or likelihood calculations should replace it in an implementation.

## What would count as two discoveries

Discovery A is a certified first-mode neighborhood. Its value is not merely estimating an eigenvector: it enables a physically executable null probe and shrinks the second-stage action set. Discovery B is a certified second-mode neighborhood obtained using that null probe; it enables a double-null probe for the third parameter. The held-out test must use a new frame from the same finite class and charge every calibration and failed comparison.

The required positive result is an amortized inequality such as

\[
N_A+N_{B|A}+N_{C|A,B}
<N_{A}+N_{B}^{\rm no\ deflation}+N_C^{\rm no\ deflation},
\]

at fixed confidence and prediction error, with the strict saving persisting on a fresh frame. The comparison must include a stateful conventional eigenmode learner allowed to perform the same adaptive deflation. If that baseline ties, the result demonstrates a useful protocol but no special advantage for a learned representation.

## Critical limitation

The proposed saving can be an artifact of restricting the baseline. A standard adaptive PCA/eigenmode learner already learns a dominant mode and deflates it. The representation test is therefore meaningful only if the learned object transfers across a new covariance instance or control family and lowers acquisition cost after its discovery cost is amortized. A one-instance eigensolver is not evidence of compounding scientific capability.

Finite candidate sets also make identifiability explicit: frames producing the same scalar-action distributions are operationally equivalent. No claim is made for arbitrary dense orthogonal frames. Correlated shots, imperfect probe normalization, and source drift must be included in the confidence budget. If the first-mode angle cannot be certified tightly enough to suppress leakage, the second discovery must abstain rather than report a false acquisition saving.
