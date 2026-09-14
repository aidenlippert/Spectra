# Approximate-commutant obstruction to a hidden density basis

This gives a Hamiltonian-only obstruction to approximating a real quartic
interaction by *any* orthogonally rotated density-density interaction. It is a
negative structural test: it cannot discover a basis, but it can rule out the
entire route below a certified tolerance.

Represent the number-conserving quartic interaction (V) as a real linear
operator on the two-particle wedge basis (ascending creation and annihilation
pairs). Fermionic sign conventions only change this representation by signed
permutations, so its Frobenius norm is invariant. Define

\[
 C_V(A)=[Q(A),V],\qquad Q(A)=\sum_{pq}A_{pq}a_p^\dagger a_q,
\]

on the real symmetric one-body matrices, with Frobenius-orthonormal basis
(E_{ii}) and ((E_{ij}+E_{ji})/\sqrt2). Let

\[
 s_1\leq\cdots\leq s_{M(M+1)/2}
\]

be the singular values of this exact linear map.

## Bound

Let (V_0) be any real orthogonally rotated density-density interaction.
Its rotated occupation projectors (P_1,\ldots,P_M) are Frobenius-orthonormal
and satisfy ([Q(P_i),V_0]=0). For (R=V-V_0),

\[
 C_V(P_i)=C_R(P_i).
\]

On a two-particle wedge basis state, an individual orbital occupation operator
has eigenvalue 0 or 1. A quartic number-conserving term changes the occupation
of at most four one-particle modes, so the sum over (i) of squared projector
eigenvalue differences is at most 4. Consequently,

\[
 \sum_{i=1}^M\|C_V(P_i)\|_F^2
 =\sum_i\|C_R(P_i)\|_F^2
 \leq 4\|R\|_F^2.
\]

Ky Fan's minimum principle gives the basis-independent obstruction

\[
 \boxed{\quad
 \|V-V_0\|_F\geq {1\over2}\sqrt{\sum_{j=1}^M s_j^2}\quad}.
\]

If (R) is represented in a coefficient basis whose individual CAR monomials
have unit Frobenius norm, then (|R|_F\leq\|R\|_1). Thus a certified lower
bound above a target (\rho) rules out every rotated density-density model
with coefficient (L^1) residual at most (\rho). One-body terms must be
handled separately: the statement applies to the quartic (V); a one-body
rotation can be transported exactly or added as a separately bounded residual.

## Exact certification route

For rational (V), use the square-root-free domain basis (E_{ii}) and
(E_{ij}+E_{ji}) and form its rational Gram (G=C_V^T C_V). The domain
Frobenius metric is (D=diag(1,2)), with 1 on diagonal generators and 2 on
off-diagonal generators. The normalized Gram is congruent to (G-tD) after
subtracting the threshold (t). To prove that the sum
of the (M) smallest squared singular values exceeds ((M-k)t), certify by
rational LDL/inertia that (G-tD) has at most (k) eigenvalues below zero.
Then

\[
 \|V-V_0\|_F\geq\tfrac12\sqrt{(M-k)t}.
\]

The inertia statement must include exact zero handling or a rational margin;
floating-point eigenvalues are a discovery diagnostic only. This test is
restricted to real orthogonal orbital rotations and the stated quartic wedge
representation until the additional block below is included. Neither version
rules out higher-body terms, number-ideal reductions, other residual norms, or
a larger certificate family. A coefficient Frobenius lower bound is not a
lower bound on ground-energy error or on the fixed-particle many-body operator
norm; it specifically rules out small coefficient residuals.

## Complex-unitary extension for real interactions

The same obstruction extends to arbitrary unitary orbital rotations when the
input quartic (V) is real in the chosen orbital basis. Use the full Hermitian
one-body domain, decomposed as

\[
\mathrm{Herm}(M)=\{A=A^T\}\oplus\{iB:B=-B^T\}.
\]

With Frobenius-normalized bases (E_{ii}),
\((E_{ij}+E_{ji})/\sqrt2\), and (i(E_{ij}-E_{ji})/\sqrt2\), the commutator
images lie respectively in real antisymmetric and imaginary symmetric wedge
matrices. These output spaces are orthogonal, so the full exact Gram is the
direct sum of the two rational Gram blocks (or their diagonal-metric versions
if square roots are cleared). The domain dimension is (M^2); select the
relevant (M) smallest singular values from the combined spectrum.

Any (M) unitary-rotated occupation projectors are Hermitian,
Frobenius-orthonormal, and commute with a unitary density-density (V_0).
The occupation-difference argument is unchanged after the complex rotation:
each quartic matrix element changes at most four occupation labels, so the
same factor-4 bound applies to the complex Frobenius norm. Therefore the boxed
bound and rational inertia certificate extend to complex unitary rotations,
provided both real and imaginary commutant blocks are included. Omitting the
antisymmetric block would incorrectly leave a complex-unitary loophole.

## Quotienting the one-body number ideal

The obstruction can also be made insensitive to a body-one number-ideal
multiplier. On the two-particle wedge, the quartic component of
\(\hat N Q(A)\) is the one-body lift
\[
L(A)=A\otimes I+I\otimes A\quad\text{restricted to }\wedge^2.
\]
Project the quartic interaction orthogonally off the range of \(L\). For a
real symmetric \(A\), diagonal generators \(L_i=L(E_{ii})\) obey
\(\langle L_i,L_i\rangle=M-1\) and \(\langle L_i,L_j\rangle=1\) for
\(i\ne j\). Off-diagonal symmetric generators have norm squared
\(2(M-2)\) and are orthogonal to the diagonal sector. If
\(t_i=\langle L_i,V\rangle\) and \(t_{ij}\) is an off-diagonal inner
product, the projection coefficients are
\[
x_i=\frac{t_i-\sum_jt_j/[2(M-1)]}{M-2},
\qquad x_{ij}=\frac{t_{ij}}{2(M-2)}.
\]
The implementation constructs the lifts by actual CAR multiplication and
checks orthogonality exactly, including every imaginary-antisymmetric lift.
For real Hermitian input their projection coefficients vanish; projecting the
real symmetric block therefore evaluates the full Hermitian projection.
The projection \(P(V)=V-L(A)\) is orthogonal and unitary
equivariant. Thus if \(V=V_0+L(X)+R\), then
\(P(V)=P(V_0)+P(R)\) and \(\|P(R)\|_F\leq\|R\|_F\). Applying the commutant
obstruction to \(P(V)\) rules out small residuals even after optimizing over
an arbitrary Hermitian body-one multiplier. The complex proof uses the full
Hermitian generator space. This section covers body-one ideals only, not
body-two or general cubic SOS multipliers.

The target remains in the density class after projection: in the density
basis of \(V_0\), every off-diagonal lift is wedge-off-diagonal and has zero
inner product with \(V_0\). The removed part is therefore a diagonal lift
\(L(\operatorname{diag}x)\), leaving \(P(V_0)\) density-diagonal. An arbitrary
off-diagonal \(X\) causes no loophole because \(P(L(X))=0\) before taking the
commutator. This argument is independent of whether \(X\) commutes with the
target occupation operators.
