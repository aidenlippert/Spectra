# Number-sector reduction of the mixed cubic dictionary

On a fixed (N)-particle sector, the mixed charge-(-1) dictionary

\[
 \{a_i\}\cup\{a_k^\dagger a_j a_i\}
\]

contains exact number-ideal dependencies. Since annihilating first leaves an
((N-1))-particle state,

\[
 \sum_j n_j a_i=(N-1)a_i,
 \qquad n_j=a_j^\dagger a_j.
\]

The (j=i) term vanishes, so every linear word can be represented by a sum of
the cubic words in the same mixed block whenever (N>1). In the charge
(+1) block the same-sector identity is
\(\sum_{j\ne i}a_i^\dagger n_j=N a_i^\dagger\), with denominator N,
not N-1. Taking an adjoint alone changes which sector is annihilated.
This is a quotient-space reduction,
not a claim that the raw words are linearly independent.

More generally, lower-degree factors can be lifted without changing the SOS
cone on the sector. If (Q) has charge zero, then

\[
 \sum_i(a_iQ)^\dagger(a_iQ)
 =Q^\dagger\Big(\sum_i n_i\Big)Q
 =NQ^\dagger Q.
\]

If (B) has charge (-2), then (a_i^\dagger B) acts on an (N-2)-particle
state after (B), and therefore

\[
 \sum_i(a_i^\dagger B)^\dagger(a_i^\dagger B)
 =B^\dagger\Big(\sum_i a_i a_i^\dagger\Big)B
 =(M-N+2)B^\dagger B.
\]

The factors (1/N) and (1/(M-N+2)) are exact rational scalings when their
denominators are nonzero. These identities show that a full mixed cubic SOS
cone contains lifted copies of the charge-zero quadratic and charge-(-2
pair) cones. Conversely, removing redundant explicit linear words from a
mixed block preserves the represented cone modulo a number-ideal multiplier;
keeping them is numerically harmless but can worsen conditioning.

More explicitly, the difference between a linear word and its cubic
representation is \(-a_i(\hat N-N)/(N-1)\) for annihilation and
\(-a_i^\dagger(\hat N-N)/N\) for creation. Substituting these relations
into a quadratic Gram form changes it by \((\hat N-N)X\), where X is
Hermitian, number conserving, and has degree at most four. Thus the reduction
requires the body-two multiplier basis. The cubic Gram is obtained by a
congruence, so it stays positive semidefinite. Conversely the reduced
dictionary is a subset of the original dictionary.

The implementation also optionally fixes the identity residual to zero. If
the equality has residual coefficient r0 on I, replacing b by b+r0 and
setting r0 to zero preserves feasibility and never worsens b-||R||1. This
removes a flat direction of the residual objective without changing its optimum.

The ordering proof is an operator identity on the (N)-sector: insert
\(\hat N-N\) times the appropriate lower-degree polynomial to rewrite the
factor expression, and replay the remainder with CAR normal ordering. It is
not an identity on the full Fock space. A verifier must retain the number
ideal and reject the reduction at (N=0) for the linear lift, or at
\(M-N+2=0\) (which is outside ordinary (0\le N\le M), but should still be
checked defensively). Boundary sectors (N=0,1) can have vanishing or
degenerate dictionaries and require explicit refusal or a separate basis.

The practical benefit is removal of (M) redundant linear directions per
charge block, while the cubic support can still grow as (O(M^3)) and the
lifted factors gain coefficients proportional to (N) or (M-N+2). Thus this
is an exact quotient and conditioning route, not evidence that the mixed SOS
cone scales polynomially for generic chemistry.
