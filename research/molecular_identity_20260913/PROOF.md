# Exact scope of the molecular local-operator pass

For distinct spin-orbital labels `p,q,r`, let `n_i = a_i† a_i`. CAR gives

\[
(a_p a_q a_r)^\dagger(a_p a_q a_r)=n_p n_q n_r\succeq0,
\qquad
(a_p^\dagger a_q a_r)^\dagger(a_p^\dagger a_q a_r)
=(1-n_p)n_qn_r\succeq0.
\]

The eight analogous products of `n_i` and `1-n_i` are mutually orthogonal
occupation projectors and sum to identity. These statements hold on the entire
fermionic Fock space, with no assumption about how many particles occupy a
local cluster. Any injective relabeling of the orbitals preserves the identities.
No determinant basis or matrix diagonalization is needed to establish them.

The previously accepted H4 restricted dual assigns

\[
y(n_2n_4n_5)
=-\frac{87413147645264941}{7000000000000000000}<0.
\]

Thus this dual cannot be moments of a physical state and is excluded by adding
the corresponding cubic square. The standard-library replay first verifies its
normalization, coefficient box, full restricted Gram PSD conditions, and all
body-two global-number ideal equalities. It then verifies 32 stored rational
negative directions from the local scan. Exact arithmetic, rather than a
floating eigenvalue, accepts each separation.

Excluding this particular dual does **not** prove that the optimized lower
bound improves: another dual feasible for the enlarged relaxation can achieve
the same objective. The three-orbital experiment illustrates why actual primal
certificate export must accompany this diagnostic.

## Coupling a short operator to an existing block

Let `B` be an existing column of same-charge, same-symmetry operators and `A`
a proposed local cubic polynomial. A positive Gram matrix on `[B; A]` permits

\[
[B;A]^\dagger
\begin{bmatrix}Q&C\\C^\dagger&D\end{bmatrix}
[B;A]
=B^\dagger QB+B^\dagger CA+A^\dagger C^\dagger B+A^\dagger DA.
\]

Keeping separate PSD blocks enforces `C=0`. The new coupled dictionary permits
these cross terms. Their support can extend beyond the four-orbital support of
`A` because the existing commutator generators can involve other orbitals. Every
cross product is normal ordered with full CAR; degree-six terms are retained.
The locality claim concerns the **added generators**, not every resulting square
or the entire Hamiltonian.

The numerical selector searches negative eigenvectors of local moment blocks.
It rounds coefficients to rationals with denominator one million. It removes
the eigenvector's linear part because all linear operators of the matching
charge/symmetry already lie in `B`: this preserves the enlarged span. It adds
the cubic part and its exact adjoint to the matching global blocks. Numerical
violations rank proposals; they are never accepting evidence for an energy.
If no existing block carries a proposed charge/symmetry, a new PSD block is
created for that channel. For example, an alpha creation together with two beta
annihilations has total charge minus one but cannot match a linear annihilator's
spin charge. Such operators still give valid number-conserving squares.

## Accepted energy interval

The unchanged symbolic verifier reconstructs

\[
H=bI+\sum_\ell F_\ell^\dagger F_\ell
+(\widehat N-N)X+R.
\]

`X` is Hermitian and number conserving. On the full fixed-total-`N` sector the
ideal term vanishes. Every canonical CAR word has operator norm at most one,
so `||R|| <= sum_w |R_w|` and

\[
E_0(H|_N)\ge L=b-\sum_w|R_w|.
\]

The Gram exporter produces rational factor rows. PSD clipping or numerical
coefficient errors can increase the exact residual and weaken `L`; they cannot
make an unverified numerical objective an accepted lower bound. The complete
rational identity and residual are rebuilt from the exported factors.

A separate existing rational state witness gives a Rayleigh upper bound `U`.
Its Hamiltonian and particle number are checked against the certificate. Hence
`L <= E_0 <= U`. Those state witnesses were discovered by FCI in finite spin
sectors; their discovery cost is explicitly retained as a reference cost.

These are energy intervals for the frozen rational STO-3G Hamiltonians. They
contain no basis-set, geometry, nuclear-motion, temperature, or experimental
error allowance. No theorem here establishes fixed-width accuracy, favorable
asymptotic cost, novelty of local positivity, or the wider mission's physical
design and preparation capabilities.
