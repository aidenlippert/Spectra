# A missing consistency condition, and a test of collective elimination

The main result is an exact positive operator that the complete old
492-generator family cannot certify as nonnegative. Its accepted moment
witness assigns that operator a negative expectation. This demonstrates
mathematical power absent from the old family. The operator is an instance
of established T1 positivity, so the contribution here is a concrete
separation and diagnosis of our restricted construction, not a new physical
law or a demonstrated scalable molecular algorithm.

The independent elimination investigation produced a compact exact response
for a known symmetric spin family and explicit failures of its simplest
transfer assumptions on H6. We kept the original fitted results and
verifier intact. No contraction fit or adaptive enrichment was run.

## The missing condition is collective three-particle positivity

For any three-form c, let

\[
C=\sum_{i<j<k}c_{ijk}a_i a_j a_k,\qquad
P=C^\dagger C+CC^\dagger\succeq0.
\]

Its expectation is the sum of two squared norms. The sixth-degree terms
cancel exactly under fermionic anticommutation, leaving only existing
one- and two-particle moments. This is the known T1 condition, explicitly
described in [Mazziotti's equations 22 and 40–42](https://arxiv.org/html/1207.0541v1)
and the earlier [three-index representability study](https://optimization-online.org/wp-content/uploads/2003/10/760.pdf).

Our rational operator has twelve terms on spatial orbitals 0, 2, 3, 5
(zero-based, interleaved spins). Writing A(i,j,k)=a_i a_j a_k:

```text
100 C = -3 A(0,1,4) - A(0,4,5) -6 A(0,6,7)
        +21 A(0,6,11) +39 A(0,7,10) -3 A(0,10,11)
        -58 A(1,6,10) +20 A(4,6,7) -5 A(4,6,11)
        -43 A(4,7,10) +13 A(4,10,11) +47 A(5,6,10).
```

The [exact separator replay](/Users/aidenlippert/Documents/Spectra/results/response_consistency_20260913/separator_receipt.json)
gives

\[
y(P)=\frac{-189044238232068984960840379}
{110000000000000000000000000000}
\approx-0.0017185839839279<0.
\]

This is a dimensionless moment, **not an energy improvement in hartrees**.
The expanded polynomial has one constant, sixteen quadratic terms and 132
quartic terms. All 149 moments it uses are explicitly present in the parent
witness; no missing values, including higher moments, were supplied.

The [fresh full-family audit](/Users/aidenlippert/Documents/Spectra/results/response_consistency_20260913/parent_family_audit.json)
accepted the unchanged parent hash
`76860d41dcdde304865ca8b7a03ca201fc3bae4afeb19f4b586c7cce88fa7206`.
It checked the quadratic rules, all 79 number-ideal equalities,
normalization, Hermiticity, the coefficient box, the fixed tail, and all
eight completed-spin blocks, of dimensions 30, 30, 93, 93, 93, 93, 30, 30.
Therefore P>=0 is not implied by their full rules. Any physical state or
admissible higher-moment extension would give y(P)>=0, which is impossible.

The support result is sharper than a negative floating eigenvalue. An
[exact support audit](/Users/aidenlippert/Documents/Spectra/results/response_consistency_20260913/support_audit.json)
proved that every T1 matrix on any three of the six original spatial
orbitals is positive definite: twenty 20x20 rational matrices. A
four-orbital relation detects the defect. The selected three-form also has
exact one-particle support rank eight, so this particular C cannot be
written using fewer than eight modes after an invertible orbital change.
This does not exclude a different separator on rotated orbitals.

The violation describes a class, not only one point. Let tau be the uniform
physical trace on the N=6 sector. The mixture (1-t)y+t tau still satisfies
every old rule by convexity, and the frozen C rejects it whenever

\[
0\le t<\frac{2392965040912265632415701}
{256392965040912265632415701}\approx0.009333193056.
\]

The predeclared 0.1% and 0.5% mixtures fail; the 1% and 2% mixtures pass
this particular condition. Three fresh rational physical trial states
also reproduced the sum-of-squared-norms identity using independent
occupation-bit action. These mixtures and physical controls establish
robustness and soundness; they are not independent molecular transfer.

## Independent molecular transfer

A separately generated H6 chain at 1.6 Angstrom was held out from selection
of the operator above. The predeclared test used exactly one solve of the
fixed full old family. It returned `user_limit`. The requested numerical
solver limit was 45 seconds, but the measured `problem.solve` call,
including its setup, took 153.695 seconds; the entire discovery subprocess
took 182.449 seconds. The exported primal lower was very weak, with a
10.785640 Ha residual penalty. That result is retained, not treated as an
energy improvement.

The resulting moment proposal, repaired with a 0.0001% physical trace
mixture, gives the frozen operator the **positive** exact expectation
2840319216396387597191 / 17187500000000000000000, approximately
0.165254936227. The predeclared numerical search also found no T1 violation
on three- or four-spatial supports; its smallest eigenvalues were 0.145269
and 0.127051 respectively. No fresh separator was exported. Full old-family
acceptance and the frozen operator's rational evaluation are recorded in
the [fresh-case replay](/Users/aidenlippert/Documents/Spectra/results/response_consistency_20260913/fresh_h6_1p6/exact_receipt.json).
That replay accepted every old-family rule in 138.755 seconds. Its dual
ceiling is above the known HF variational upper, so this fresh run did not
reproduce the original energy obstruction either.

**Independent molecular transfer remains unestablished.** This bounded,
poorly converged solve did not produce a new T1-rejected counterexample.
That is not evidence against the universal identity, nor a proof that no
other fresh witness or larger-support constraint separates. The numerical
support scan is not promoted to an exact all-support positivity theorem.
No solver retries, T1 energy optimization, or FCI calculation were used.

## A response that eliminates a large sector, and where it stops working

Starting with the supplied block congruence, choose X and compute
E=B*−DX and K=A−BX−X*B*+X*DX. The sufficient conditions
D>=delta I>0 and K>=E*E/delta are an exact completion of squares. This
fits established [Feshbach-Schur methods](https://arxiv.org/abs/2105.02058).
The full derivations, including a conditional charge-shift extension,
are in [DERIVATION.md](/Users/aidenlippert/Documents/Spectra/research/response_consistency_20260913/DERIVATION.md).

For a central spin uniformly coupled to M bath spins at fixed total
excitation q, take

\[
H=\epsilon n_c+\omega N_b+\chi n_cN_b
 +g(s_c^+J_-+s_c^-J_+).
\]

The two diagonal blocks are scalars A0=omega q and
D0=epsilon+(omega+chi)(q−1). One collective response,

\[
X=\frac{g}{D_0-b}J_-,
\]

has E=0. The exact identity

\[
q(M-q+1)I-J_+J_- = \sum_{i<j}(I-\mathrm{Swap}_{ij})\succeq0
\]

reduces the lower certificate to
(A0−b)(D0−b)>=g²q(M−q+1), with D0−b>0.
A symmetric Dicke trial gives a rational upper without listing its
amplitudes. The scalar gap, response and remainder certificates are all
explicit. This realizes collective elimination even when coupling is
comparable to detuning; it requires exact homogeneity. Homogeneous
central-spin solutions are [established mathematics](https://arxiv.org/html/1810.03012v4).

The [exact response replay](/Users/aidenlippert/Documents/Spectra/results/response_consistency_20260913/response_receipt.json)
tested nine formula-specified cases, from two to one million bath spins.
At one million spins the model plus certificate occupies 270 JSON bytes
and the exact interval width is below 5.59e−10 in model energy units.
No many-body matrix, inverse, amplitudes, or binomial coefficient was
constructed on this accepting path. The input's homogeneous formula is
essential: recognizing that structure in a supplied list of couplings
still costs at least the input length. Applying a compact operator formula
to an arbitrary explicit many-body vector is not thereby cheap.

The [H6 falsification](/Users/aidenlippert/Documents/Spectra/results/response_consistency_20260913/h6_closure.json)
found all 45 pattern pairs noncommuting. Their exact commutator closure
has dimension 35, the entire traceless 6x6 matrix algebra. For the concrete
split by occupation of spin orbital zero, two determinant diagonals in
the eliminated block differ by exactly 0.180511162346 Ha; D is not scalar.
These refute the tested commuting/scalar closure assumptions. A
35-dimensional one-particle algebra is still polynomial in size, so this
is not an impossibility theorem for compact noncommutative response.

An [independent inhomogeneous spin control](/Users/aidenlippert/Documents/Spectra/results/response_consistency_20260913/response_controls.json)
also exposes the boundary. Changing one of six couplings from 1/6 to 1/4
gives an exact trial upper −0.4930842352913686, below the homogeneous
lower −0.4465176722724529. Carrying over that bound would be false.
The checker refuses the changed model. The scalar-D elimination still
exists for weighted couplings; it is the old collective norm certificate
that no longer applies. Other extensions are not ruled out.

## Evidence, cost and remaining scientific question

Original separator discovery built a 220x220 T1 moment matrix and took
2.886 seconds. A 3.015-second reproducibility run records every rounding
attempt, including the unsuccessful denominator-10 proposal; denominator
100 reproduces the same separator. This search has not demonstrated a
scalable selection algorithm. The fresh parent audit took 180.347 seconds,
and exact three-orbital/support validation took 0.616 seconds. Additional
stage timings and the independent molecular test are in
[cost_ledger.json](/Users/aidenlippert/Documents/Spectra/results/response_consistency_20260913/cost_ledger.json).
These are measured computation times, not total research/calendar time.
Including the unsuccessful transfer experiment and both test runs, the
sum of measured stage times is 510.210 seconds. The small exact stages
overlapped, so this sum is not a claim about calendar elapsed time.

The new spin controls enumerated 209 basis states in small diagnostic
matrices. The independent fermionic controls generated 1,914 basis labels
to select sparse trial vectors; H6's nonscalar-D test evaluated two
determinants. These costs are separate from the nonenumerating accepting
rules. Fresh molecular fixture generation uses an HF control only. All
numerical proposals remain untrusted until rational replay.

The six focused tests cover exact solvable cases, gap/remainder refusal,
refusal of a different Hamiltonian, witness hash binding, rejection of
higher moments, and the elementary occupancy identity. All passed.
Accepting replays ran under Python `-S` and confirmed no NumPy, SciPy,
CVXPY or PySCF imports. All 665 files in the preceding frozen manifest
remain unchanged, including the old verifier.

This pass establishes an explicit missing consistency condition and
identifies the assumptions behind one compact elimination
formula. It does not establish that a few patterns capture all molecular
correlation, or that the added condition improves a molecular energy bound.
The old witness's 1.513465 mHa fixed-upper floor still describes the old
family only; removing that witness does not calculate a new optimum.

The next foundational target is a useful nonscalar response with a cheap
independent gap certificate and a collectively controlled remainder.
The charge-shift identity in the derivation specifies one sufficient
structure to test. The verified T1 defect is now a concrete consistency
requirement for proposed compressions. Further fitting inside the old
span cannot supply that missing condition.
