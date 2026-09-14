# Dual obstruction plan for the commutator frame

The enriched H-only dictionary adds charge-(-1) factors
\([H,a_i]\), whose normal-ordered content is cubic plus lower degree. The
observed H4 lower values are indistinguishable from the quadratic baseline,
but this is not an algebraic redundancy theorem. In general
\([H,a_i]\) is not in the quadratic factor span modulo the body-two number
ideal, and it can enlarge the primal SOS cone. A zero numerical gain is
therefore evidence about pricing/conditioning, not proof of equality of cones.

## Exact dual certificate of a missed lower bound

Write the restricted primal in canonical number-conserving word coordinates
as

\[
 f=b\,1+\sum_r x_r I_r+\sum_\alpha G_\alpha(M_\alpha),
 \qquad G_\alpha\succeq0,
\]

where (I_r) are the body-two number-ideal rows and (M_\alpha) is the Gram
polynomial of each retained factor block. To prove that a target (b_*>b) is
outside this cone, construct an exact dual functional (y) on the full
canonical coefficient rows satisfying:

1. (y(1)=1) and (y(I_r)=0) for every ideal row;
2. (y(f_H)=b_y), with (b_y<b_*) (or normalize the objective row directly);
3. (y(M_\alpha)\succeq0) for every restricted Gram block;
4. (y) respects all CAR canonicalization and Hermitian pairing constraints.

Then every feasible restricted certificate has (b\le y(H)), because ideal
terms vanish and each PSD Gram contributes a nonnegative pairing. A strict
gap (y(H)<b_*) proves the restricted cone cannot reach (b_*), independent
of floating-point optimizer failure.

The implementation form is a rational feasibility problem. Assemble the
linear coefficient map (T) from ideal and Gram variables to canonical rows,
solve the free-column constraints (F_free^T y=e_b) together with (|y_w|\le1), and
certify each dual Gram (\sum_w y_w\,\mathrm{GramPoly}_\alpha(w)\succeq0) by
exact rational LDL/congruence. For the implemented coefficient-L1 residual objective, the box constraint
|y_w|<=1 is mandatory. Omitting it would only obstruct exact zero-residual
identities, not the certificates actually optimized here.

## Maximally mixed fixed-sector seed

A useful dual seed is the maximally mixed (N)-particle state on (M) modes.
For a normal-ordered balanced word whose creation and annihilation supports
are the same (k)-set, its moment is

\[
 \omega_k=\frac{\binom{N}{k}}{\binom{M}{k}},
\]

times the CAR sign from canonical ordering; words with unequal supports have
zero moment. This functional is positive on every square (A^\dagger A) and
has exact rational entries. It already satisfies the number-ideal affine constraints exactly. A rounded
optimizer candidate must first be corrected onto those constraints; a small
convex mixture with this exact seed can then restore positivity on the
quotient by exact Gram null directions. Such a repair must be proved feasible,
not inferred from floating eigenvalues.

For the H4 restricted experiment (M=8,N=4, 777 canonical rows and a 111-row
dual basis), the concrete audit should export (y), all ideal residuals,
each dual Gram's exact LDL inertia, and the resulting (y(H)). A successful
strict inequality against the enriched primal lower target is a proof that the
restricted cone misses that target. Failure of this construction leaves the
question open; it does not establish redundancy. The same recipe applies to
the H6 enriched run, with the larger rational row map and exact bit-size
controls.

The comparison should also retain the independently certified physical upper
bound. A dual obstruction only separates a restricted lower cone from a target
certificate; it says nothing by itself about the true ground-state energy or
about the usefulness of the commutator frame for other Hamiltonians.

For a physical accuracy obstruction, compare the dual objective with a
rigorous lower bound on the true ground energy, such as the full H4 cubic
certificate lower -3.667000108966591573. If that lower minus y(H) exceeds
.0016, no certificate in this restricted L1-penalized cone can have a true
ground-energy error <=.0016. A physical upper witness alone is insufficient
for this stronger statement. No such dual witness has yet been constructed.
