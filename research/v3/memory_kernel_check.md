# Finite dissipative block certificate

The executable model is explicitly

\[
L=\begin{pmatrix}K-D_P&B\\-B^T&-D_Q\end{pmatrix},
\]

with real skew `K`, positive diagonal `D_P,D_Q`, and zero initial discarded coordinates. Let `gamma` be the smallest rate in **both** sectors, let `b >= ||B||_2`, and let `q >= ||x(0)||_2`. The symmetric part of the full generator and each principal block is bounded above by `-gamma I`. Therefore all their propagators contract by `exp(-gamma t)`.

Variation of constants gives

\[
\|y(t)\|_2\le b q t e^{-\gamma t},
\qquad
\|z(t)-\hat z(t)\|_2\le\frac{b^2q t^2}{2}e^{-\gamma t}.
\]

The stored finite-time bound drops the exponential factor; the stored uniform bound uses the conservative inequality `t^2 exp(-gamma t)/2 <= 2/gamma^2`. Rational square-root upper bounds give `b <= sqrt(||B||_1 ||B||_infinity)` rounded upward, and bound the initial norm. The checker verifies all dimensions, rates, skew symmetry and arithmetic inequalities.

Damping only the discarded block does not establish this uniform retained-error theorem. That omission was found during review and corrected: retained damping is now explicit and zero rates are rejected. The result also excludes nonzero initial discarded coordinates, which would add an inhomogeneous transient term.

These are retained coefficient 2-norm bounds. They do not by themselves bound an arbitrary physical observable without a dimension-dependent norm conversion. The separate [local Pauli coefficient theorem](local_pauli_theorem.md) supplies a physical-error route under a stronger damping assumption.
