# A retained interference system with a collectively certified bath

This proof applies to the declared model family. It does not establish that a
molecular Coulomb Hamiltonian has this representation. The accepting program
is `collective.py`; numerical diagonalization only proposes its witnesses.

## 1. Hamiltonian and fixed-number identity

There are four active fermion modes, and two bath bands with B = L² modes each.
Write N_A and N_B for active and total bath particle number. Work in the sector
N_A + N_B = B + n₀, with n₀ = 2 in the scaling campaign. Let

\[
e_j=g+w\frac{j}{B-1},\qquad j=0,\ldots,B-1,\qquad g>0,\quad w\ge0,
\]

with e₀ = g when B = 1. Positive/negative bath one-body energies are ±e_j.
The active Hamiltonian H_A is arbitrary within the checker's four-mode,
real, number-conserving, degree-at-most-four input grammar. The campaign uses
nonzero hopping on triangles, including a frustrated sign product, and
nonzero two-body density interactions; its coefficients are in each model.

For the two specified active combinations A₊ and A₋, the hybridization is

\[
V=\frac{t}{L}\sum_{\sigma\in\{-,+\}}
\left(A_\sigma^\dagger\sum_{j=0}^{B-1}b_{\sigma j}
+\sum_{j=0}^{B-1}b_{\sigma j}^\dagger A_\sigma\right).
\]

The 1/L normalization is part of the input family. It holds the total
hybridization norm bounded as B grows. Fixed per-edge coupling would have a
different size dependence and is not covered by the reported scaling result.

The model also contains extensive collective charging interactions

\[
u(N_A-n_0)^2+v(N_B-B)^2+z(N_A-n_0)(N_B-B).
\]

On the specified total-number sector these equal
(u+v−z)(N_A−n₀)². Define H̄_A by adding this term to H_A. This is an operator
identity on the sector: it does **not** assume that N_A has a definite value in
an interacting eigenstate. The Hamiltonian conserves total particle number;
hybridization still changes N_A.

Set E_fill = −Σ_j e_j, calculated exactly by the arithmetic-series formula.
Then, on this sector,

\[
H=E_{\rm fill}+\bar H_A+V+
\sum_j e_j n_{+j}+\sum_j e_j(1-n_{-j}). \tag{1}
\]

The filled negative bath is a reference for rewriting the operator, not an
assumption about actual occupations.

## 2. Operator envelopes keep all coupling

Partition each band into bins. In each bin I let e_I⁻ = min e_j and
e_I⁺ = max e_j. Form H⁻ and H⁺ by replacing the coefficients of the particle
and hole terms in (1) by the corresponding endpoints. Keep E_fill, H̄_A, and
V identical in both envelopes. Because n_j and 1−n_j are positive semidefinite,

\[
H^-\preceq H\preceq H^+. \tag{2}
\]

This order holds irrespective of commutators between the common interactions,
couplings, and bath terms. Replacing the extensive filled-bath energy by an
endpoint approximation would lose the useful size dependence; (2) preserves
its exact value. Refining a partition improves the **optimal** envelope
endpoints by the same operator order. Finite-precision proposals are still
subject to their separate exact checks.

## 3. Exact elimination inside each envelope

Within a bin of γ modes the envelope bath energies are equal. Its normalized
bright mode is β_I = γ⁻¹ᐟ² Σ_{j∈I} b_j. Complete this to a unitary basis of that
bin. Only β_I couples to the active modes. The other γ−1 modes are dark and
appear only in the envelope's diagonal bath excitation energy. Collective
charging terms have already been reduced by the fixed-number identity.

The code retains k = 4 + 2q modes for q bins per band. It does not build the
unitary matrix or enumerate dark modes. For each retained particle number n,
the dark subsystem must contain K = B+n₀−n particles. If K is outside
[0, Σ_I(γ_I−1)], the sector is infeasible. Otherwise the dark excitation energy
is the minimum of

\[
\sum_{I\subset -} e_I(\gamma_I-1)
+\sum_I\sigma_I e_I k_I,
\quad \sum_I k_I=K,\quad 0\le k_I\le\gamma_I-1. \tag{3}
\]

All k_I are integers. Filling the sorted one-body costs σ_I e_I by their
multiplicities gives this minimum exactly. At most 2q costs are sorted.
For fixed n, the same dark minimum applies to every retained vector because
the retained Hamiltonian preserves n. No omitted dark state can lower (3).

To avoid irrational square roots in acceptance, use unnormalized bright
operators c_I = Σ_{j∈I} b_j, whose CAR are {c_I,c_J†} = γ_I δ_IJ.
For the occupation label s of a retained basis vector, its norm is

\[
G_s=\prod_{I\ {m occupied\ in}\ s}\gamma_I.
\]

Creation in this basis has the ordinary fermionic sign; annihilation also
has the factor γ_I. The form matrix is M = G times the action matrix, so it
is exactly symmetric. The physical orthonormal matrix would be
G⁻¹ᐟ² M G⁻¹ᐟ². The checker constructs M and G directly with rational arithmetic,
including (3); it never substitutes an ordinary CAR metric for the bright
operators. These are actual many-body retained matrices, not a lifting of
one-particle positivity to an interacting occupation problem.

## 4. Accepting an energy interval

Subtract E_fill and reconstruct every feasible retained sector of H⁻ and H⁺.
For a proposed rational lower endpoint b, check

\[
M_n^- - bG_n\succeq0\quad\text{for every feasible }n. \tag{4}
\]

Exact symmetric elimination verifies this condition, rejecting negative pivots
and a nonzero off-diagonal row at a zero pivot. Since G_n is positive definite,
(4), the dark minimum, and (2) prove E₀(H) ≥ E_fill + b.

A nonzero integer retained vector x, together with its feasible particle
number and the minimizing dark occupations, gives the exact upper Rayleigh
value u⁺ = xᵀM⁺x/(xᵀGx). Thus E₀(H) ≤ E_fill + u⁺. The physical state
represented by this witness need not be written as a large determinant list.
The model hash, all input domains, particle numbers, amplitudes, metrics, and
forms are checked or reconstructed on replay. Saved floating eigenvalues or
saved endpoint receipts cannot accept a certificate.

## 5. A size-independent error ceiling

Let X count all bath particles and holes relative to the filled reference.
Writing H⁻−E_fill = H̄_A+V+D⁻ gives D⁻ ≥ g_min X. Choose certified constants
H̄_A ≥ aI and ||V|| ≤ Γ. The implementation uses the conservative bounds

\[
a=-\sum_{\text{normal-ordered terms of }\bar H_A}|h_w|,
\qquad
\Gamma=2|t|\sum_{\sigma,i}|(A_\sigma)_i|.
\]

Every normalized CAR monomial has norm at most one, and
||Σ_j b_j|| = √B = L; these inequalities justify both constants. For a
normalized lower-envelope trial state with exact shifted energy u⁻,

\[
\langle X\rangle\le\frac{u^- -a+\Gamma}{g_{\min}}.
\]

If δ is the maximum bin width, 0 ≤ H−H⁻ ≤ δX. Applying the variational
principle to this same trial state yields a second upper endpoint and hence

\[
0\le E_0(H)-(E_{\rm fill}+b)
\le (u^- -b)+\frac{\delta}{g_{\min}}(u^- -a+\Gamma). \tag{5}
\]

The checker uses the smaller of this upper endpoint and the H⁺ Rayleigh
endpoint. Equation (5) is independent of bath multiplicity at fixed trial
quality, gap, bandwidth, active Hamiltonian, and normalized coupling. It is
deliberately coarse: the default one-bin ceiling is about 2.29046 model units,
whereas the directly certified interval is about 0.01509. A small gap worsens
the ceiling. This is not a gap-independent arbitrary-accuracy theorem.

There is a stronger structural observation for the one-bin campaign. For
B ≥ 4 and n₀ = 2, every n = 0,…,6 is feasible, and (3) is e|3−n| in either
envelope. The normalized bright couplings are (t/L)√B = t. Consequently the
two normalized retained Hamiltonian families are exactly independent of B;
only E_fill and the rational basis metric change. This explains the constant
observed widths for B = 4,…,4096 without extrapolating from timings alone.

## 6. Resources and the failure boundary

For fixed q, the total number of retained basis states is at most 2^(4+2q).
The largest sector has dimension at most C(4+2q, floor((4+2q)/2)). Acceptance
uses rational matrices of those dimensions and O(q) spectator arithmetic per
sector. With one/two bins the largest matrices are 20/70, and the total
allocated entries per envelope are at most 924/12,870. Exact arithmetic and
witness bit lengths grow with log B. The implementation enforces q ≤ 2.
Increasing the number of retained modes still has exponential cost; no such
cost is claimed to disappear.

The compact model is a formula for its energies and coupling coefficients.
Reading an arbitrary explicitly listed molecular tensor would incur its
input-processing cost. The resource statement does not hide that cost or
assert that arbitrary physical inputs admit the formula.

Rank-one coupling alone is insufficient for exact dark elimination. For
distinct bath energies e₁,…,e_B and a coupling vector with all entries nonzero,
the Krylov vectors v, Ev,…,E^(B−1)v have determinant
(Π_j v_j) Π_{i<j}(e_j−e_i), which is nonzero. The smallest bath subspace invariant
under E and containing v can therefore be the entire bath. Our dispersed
construction uses (2) before invoking degeneracy; it does not discard this
bright–dark mixing from the original Hamiltonian.

## 7. Meaning of the molecular audit

For each four-orbital subset in the frozen H4/H6 orbital basis, retain all
quadratic terms and all active-only quartic terms. Fit two collective density
terms J_B N_B(N_B−1)/2 and J_AB N_A N_B by exact coefficient medians. The
remaining Hermitian quartic operator is W. The median minimizes the density
coefficient L1 cost for each disjoint class; the scan chooses the smallest
paired-monomial sufficient norm bound over those subsets.

For off-diagonal number-conserving monomials, a conjugate pair has norm at
most its coefficient magnitude; each diagonal monomial is a signed
occupation projector. Summing these costs gives ||W|| ≤ η. If an interval
for the retained model were available, an unconditional norm transfer would
widen it by 2η. The audit does not solve that retained molecular model, whose
quadratic terms need not satisfy our two-band assumptions. These budgets are
not newly computed molecular energy intervals.

For the selected W, the audit also provides two normalized occupation states
in the physical fixed-N sector with an exact nonzero matrix element c.
Then ||W|| ≥ |c|, even after an arbitrary scalar shift. Discovery scans the
small molecular sector and is charged as such; replay applies W to one state.
This lower bound concerns that particular truncation, not another choice of
patterns or a more useful order-based certificate.

Finally, a nonzero minor certifies the rank of a one-leg flattening of the
full-Fock quartic interaction tensor. Orbital rotations preserve this rank;
an interaction supported on k orbitals has rank at most k. This only obstructs
an orbital-rotation-only strict impurity representation of the **same full
tensor**. It does not include fixed-number identities or collective
interacting spectators. The control model in this pass itself has full
one-leg rank and still permits the elimination proved above. Neither this
rank calculation nor the residual budget is a no-go theorem for the original
many-electron compression hypothesis.
