# Trace-seed Chebyshev extension beyond reflection positivity

There is a valid generalization that does not require the Hubbard
particle-hole positivity theorem.  Let `H` be any Hermitian Hamiltonian on a
declared finite sector of dimension `D`, and define left multiplication on
Hilbert-Schmidt operator space by

\[
 \mathcal L(X)=HX.
\]

With inner product `\langle X,Y\rangle=\operatorname{Tr}(X^*Y)`,
`\mathcal L` is self-adjoint because
`Tr(X^*HY)=Tr((HX)^*Y)`.  Its spectrum is exactly the spectrum of `H`, each
eigenvalue repeated `D` times.  In particular

\[
 W=|g\rangle\langle g|,\qquad \mathcal L(W)=E_0W,
\]

is a PSD ground eigenmatrix for an arbitrary ground state `|g>`, and the
identity seed satisfies

\[
 \left\langle W,\frac{I}{\sqrt D}\right\rangle
 =\frac{1}{\sqrt D}.
\]

Therefore the inexact Chebyshev exclusion argument applies verbatim to
`\mathcal L`, with the seed overlap `1/\sqrt D`; no real-amplitude or
reflection-positivity assumption is needed.  Recurrence vectors need not stay
Hermitian, since only their Hilbert-Schmidt projection onto `W` is used.

## Scope and hidden cost

The result certifies the lowest energy only in the declared finite particle,
spin, symmetry sector.  It does not compare different particle numbers or
external physical-model errors.  The identity should mean the identity on
that sector, not the full Fock-space identity.

This extension is mathematically broader but not automatically computationally
efficient.  A sector-restricted identity is itself a number-conserving MPO and
may have bond dimension growing with system size.  Left multiplication by a
local Hamiltonian acts on operator space with local dimension `d^2`; an MPS
state bond `\chi` can become an operator/MPO bond of order `\chi^2` after
purification or vectorization, while a sum of Hamiltonian terms can add further
bond factors.  Every compressed recurrence must therefore provide an explicit
HS residual certificate.

The cheap acceptance route is an MPO/MPS representation of each recurrence
vector, direct residual construction, and an outward-rounded mixed-transfer
contraction of its HS norm.  Post-hoc truncation without a certified residual
does not prove the Chebyshev inequality.  A claimed scaling theorem must bound
sector-identity bond, operator bond, recurrence count, contraction cost, and
the `D^{-1/2}` seed-overlap penalty.

Thus the trace seed removes a symmetry restriction and gives a universal
finite-sector lower-bound framework.  It does not by itself solve the dense
operator-space representation problem or establish a general efficient
many-body algorithm.
