# History closure and companion equivalence for the v6 benchmark

Consider the noiseless LTI system

\[
x_{t+1}=Ax_t+Bu_t,\qquad y_t=Cx_t,
\]

with state dimension (n), input \(u_t\in\mathbb R^m\), output \(y_t\in\mathbb R^q\), and no direct feedthrough. Let

\[
\chi_A(\lambda)=\lambda^n+c_{n-1}\lambda^{n-1}+\cdots+c_0.
\]

Cayley–Hamilton gives (A^n=-\sum_{j=0}^{n-1}c_jA^j). Unrolling the state and collecting the input terms gives

\[
y_{t+n}+\sum_{j=0}^{n-1}c_jy_{t+j}=\sum_{r=0}^{n-1}D_ru_{t+r},
\]

where

\[
D_r=C\left[A^{n-1-r}+\sum_{j=r+1}^{n-1}c_jA^{j-1-r}\right]B.
\]

The derivation must collect the input contributions before identifying state terms: (CA^jx_t=y_{t+j}) is true only in the zero-input case. Reindexing gives an ARX law involving (y_t,\ldots,y_{t-n+1}) and the current and preceding (n-1) inputs. Order (n) is an upper bound; pole–zero cancellations can reduce the minimal order.

This identity is deterministic. With process or observation noise, it produces a stochastic residual containing filtered process and measurement noise. The residual can be colored and correlated with regressors, so finite-sample OLS is not certified by the noiseless identity. Statistical guarantees need separate excitation, stability, noise, and errors-in-variables assumptions.

## Exact companion equivalence

Fix (p\ge1) and arbitrary learned coefficients in

\[
y_{t+1}=\sum_{i=0}^{p-1}a_i y_{t-i}+\sum_{i=0}^{p-1}b_i^{\mathsf T}u_{t-i}+d.
\]

Define

\[
z_t=[y_t,y_{t-1},\ldots,y_{t-p+1},u_{t-1},\ldots,u_{t-p+1},1].
\]

The current input (u_t) is exogenous. Let (w=[a_0,\ldots,a_{p-1},b_1^{\mathsf T},\ldots,b_{p-1}^{\mathsf T},d]). A companion implementation sets

\[
\widehat y_{t+1}=w^{\mathsf T}z_t+b_0^{\mathsf T}u_t,
\]

then shifts the output register by inserting \(\widehat y_{t+1}\), shifts the stored input register by inserting (u_t), and keeps the constant coordinate equal to one. An order-(p) ARX model with the same coefficients performs exactly the same computation.

The proof is induction on (t). Initially both methods receive the same output/input history. Their first forecasts are identical by the displayed equation. If forecasts and registers agree through time (t), both receive the same next exogenous input, apply the same recurrence, and perform the same shifts. They therefore agree at (t+1), and hence for every future input sequence.

This holds for every coefficient tuple, including coefficients not arising from a physical realization. Thus a privileged controlled delay-state/compiler representation has no additional predictive expressivity over an order-complete ARX family with the same history and current input. It may still offer a smaller implementation, better inductive bias, easier certification, or lower sample complexity; those require separate resource or statistical evidence. Giving the privileged method the correct order or grammar while denying equivalent order selection or feature construction to ARX would create artificial headroom.

The theorem does not claim that capped, sparse, retrieval, or otherwise weaker estimators match the companion model. It removes only an expressivity-based advantage against the order-complete comparator. The v6 benchmark must therefore compare prediction, coverage, control, sample efficiency, and total computation under matched information and costs, while keeping noiseless closure separate from noisy identification.
