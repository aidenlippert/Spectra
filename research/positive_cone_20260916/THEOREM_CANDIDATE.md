# Auxiliary semigroup argument

The accepted elementary proof and complete model binding are in
[DERIVATION.md](DERIVATION.md). This note records an alternative explanation
of the positive ground eigenmatrix, not an additional computational result.

For ladder labels `i=2r+s`, the correct bipartite gauge is

\[
\eta_{2r+s}=(-1)^{r+s}.
\]

After down-spin particle-hole conjugation, in units of `t`,

\[
\mathcal L(C)=KC+CK-8\sum_i N_i C N_i+32C.
\]

The semigroup `exp(-tau L)` preserves the cone of complex Hermitian PSD
matrices. The hopping part acts by

\[
C\longmapsto e^{-\tau K}Ce^{-\tau K},
\]

which is completely positive. For `P_i(C)=N_i C N_i`, the projector identity
`N_i^2=N_i` gives

\[
e^{8\tau P_i}(C)=C+(e^{8\tau}-1)N_i C N_i,
\]

also completely positive for `tau>=0`. The constant supplies the positive
scalar factor `exp(-32 tau)`. Lie-Trotter convergence proves that the full
semigroup is positive. It need not preserve trace.

Finite-dimensional cone Perron-Frobenius supplies a nonzero PSD eigenmatrix
at the spectral radius of this semigroup. Because `L` is self-adjoint, that
radius is `exp(-tau e_0)` and the associated eigenmatrix belongs to the
lowest eigenspace of `L`. The direct absolute-value proof in the derivation
establishes the same fact without invoking this theorem.

Consequently, `C>0` and `L(C)-ell C>=0` imply `ell<=e_0` by the trace
argument. The actual rational candidate passed exact replay with width
`4.774093164360818e-7 t`.

The cone variable has dimension `70 x 70`, containing all 4,900 amplitudes.
An ordinary rank-deficient matrix does not meet the strict-positivity gate.
A compact tensor representation could still describe a full-rank matrix,
but a compact construction and exact residual-positivity check have not
been established here.
