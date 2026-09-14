# Warm control-error theorem: exact audit

Let (Q) be a fixed orthogonal (4\times4) matrix and (r_0=Qe_1). Define

\[
u_i=\frac35r_{i-1}+a_i\frac45Qe_{i+1},\qquad
r_i=-a_i\frac45r_{i-1}+\frac35Qe_{i+1},
\]

with (a_i\in\{\pm1\}). Since (r_{i-1}\perp Qe_{i+1}), each pair ((u_i,r_i)) is orthonormal and the three (u_i) are mutually orthogonal by induction. The source at stage (j) is

\[
C_j=I+\sum_{i\le j}\lambda_i u_i u_i^T,
\qquad (\lambda_1,\lambda_2,\lambda_3)=(10^{14},10^{10},10^6).
\]

Assuming the previous signs are known, the proposed (+) null probe is

\[
x_j=-\frac45r_{j-1}+\frac35Qe_{j+1}.
\]

For the correct (a_j=+1), (x_j=r_j), so it is orthogonal to all active modes and has variance exactly (x_j^TC_jx_j=1). For (a_j=-1), (x_j\cdot u_j=-24/25), giving variance (1+(24/25)^2\lambda_j). These are exact fractions.

## Actuator perturbation

Suppose the applied unit vector \(\hat x_j\) satisfies \(\|\hat x_j-x_j\|_2\le\eta=10^{-9}), independently of the hidden sign. For any unit mode (u), 

\[
|\hat x_j\cdot u|\le\eta
\]

in the null case, hence

\[
\operatorname{Var}(Y_j)\le1+\eta^2\sum_i\lambda_i
 =1+10^{-18}(10^{14}+10^{10}+10^6)<1001/1000.
\]

In the alternative case,

\[
|\hat x_j\cdot u_j|\ge24/25-eta>19/20,
\]

so the mode contribution is at least ((19/20)^2\lambda_j).

## Executed classifier and error accounting

Use the fixed-scale classifier |Y|>20 with (B=I), and let the readout add bounded error |\omega|\le1/1000. For source couplings scaled by an unknown α∈[1,2], the null variance remains below (1001/1000) (the larger α only affects active-mode terms, whose null projections are bounded by η). A null false positive implies |Y|≥19.999, so Chebyshev gives

\[
P(|Y|>19.999)\le\frac{1001/1000}{19.999^2}<\frac1{399}.
\]

For the alternative, a miss implies |Y|≤20.001. The centered Gaussian density is bounded by (1/(\sqrt{2\pi}\sigma)), so the interval ([-20.001,20.001]) of width 40.002 has probability at most (40.002/(\sqrt{2\pi}\sigma)<20.001/(.95\sqrt{\lambda_j})). This is bounded by (1/399) for (j=1,2) and (11/500) for (j=3), using λ=(10^14,10^10,10^6) and the lower projection (0.95).

Sequential conditional union accounting gives

\[
\frac2{399}+\frac{11}{500}=\frac{5389}{199500}<0.05.
\]

The union remains valid only when each stage's actuator error and measurement noise are conditionally independent of the hidden current sign, and when a wrong earlier sign is included in the conditional failure event. If Gaussian preparation has total-variation error at most (10^{-9}) per stage, add (3\times10^{-9}) to the final risk. Held-out couplings may change the mode hierarchy; the construction remains valid only if the coupling preserves the declared orthogonal frame and lower bound on the weakest active eigenvalue.
