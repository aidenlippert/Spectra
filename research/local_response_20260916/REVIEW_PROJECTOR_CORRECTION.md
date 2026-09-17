# Review of the grouped projector correction

The correction is mathematically valid under the stated rank-one and support
assumptions.

For each patch, `profile_check` certifies

\[
 A_i:=H_i-s_iI\succeq g_i(I-P_i),
\]

where `P_i` is the rank-one projector onto the submitted normalized local
vector.  The matrix assembled by the checker is exactly
`A_i-g_i(I-P_i)`, including the normalization by the integer-vector norm.

The two outer patches have disjoint supports, so their rank-one projectors
commute and `P_0P_2` is rank one globally. The middle projector has rank one
on BC and is tensored with the identity on A,D; it is not rank one globally. The
computed `p` is exactly

\[
 p=\|(P_0P_2)P_1\|^2
  =\langle u_0\otimes u_2,P_1(u_0\otimes u_2)\rangle/(n_0n_2),
\]

computed from the reduced outer densities and the middle vector.  The local
fermion sign convention does not change this scalar: each submitted vector
has definite spin populations and therefore definite fermion parity, making
the disjoint even local projectors commute under the Jordan-Wigner embedding.

The grouped operator obeys

\[
 A_0+A_1+A_2\succeq g(I-P_0P_2)+h(I-P_1),
 \quad g=\min(g_0,g_2),\ h=g_1.
\]

For arbitrary-rank projectors with operator-norm overlap squared at most `p`,
the smallest eigenvalue of the right side is bounded below by the smaller root of

\[
 (g-\gamma)(h-\gamma)=ghp.
\]

Thus the checker condition
`0<=gamma<=min(g,h)` and
`(g-gamma)(h-gamma)>=gh*p` correctly certifies the correction.  The full
lower bound is `sum(s_i)+gamma`.

The reported comparison is also logically sound: if the old-family dual
ceiling is obtained from a feasible compatible local marginal assignment,
then every old-family lower shift sum is at most that ceiling.  Since the new
accepted lower bound is at most the true ground energy, a positive difference
`new_lower-old_family_ceiling` proves that the old declared family cannot
reach the true energy below that separation.  It is a bound relative to the
fixed Hamiltonian and fixed physical upper convention; it is not itself a
physical upper bound or a universal limitation.

Required scope caveats: the result relies on rank-one local profiles, exact
disjoint-support tensor factorization (including parity), and exact local PSD
replay.  It does not automatically extend to overlapping outer projectors,
higher-rank profiles, odd local operators, or recursive composition.  Numerical
eigenvectors and proposed gaps are discovery inputs only; the exact replay is
the evidence.
