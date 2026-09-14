# Residual domination: a bounded local improvement

This experiment isolates one concrete compression rule for an omitted
number-conserving one-body CAR residual

\[
 R=d\Gamma(A)=\sum_{i<j} c_{ij}(a_i^\dagger a_j+a_j^\dagger a_i).
\]

The usual coefficient allowance charges each Hermitian word separately,
giving

\[
 \eta_{\ell_1}=2\sum_{i<j}|c_{ij}|.
\]

Let \(r_i=\sum_j|c_{ij}|\).  On the fixed \(N\)-particle sector,
\(\|d\Gamma(A)\|\le N\|A\|_\infty=N\max_i r_i\).  This follows by
writing \(d\Gamma(A)=\sum_{k=1}^N A^{(k)}\) on the antisymmetric tensor
power and applying the triangle inequality.  It is therefore a valid local
rational residual budget:

\[
 \eta_{\rm row}=N\max_i r_i.
\]

The verifier only accumulates rational row sums.  It does not construct the
global fixed-sector Fock matrix, and its arithmetic cost is `O(|E|+M)` for
`|E|` residual edges on `M` modes.  A four-mode matrix audit is included in
the script solely to independently check the inequality on the fixture.

## Result

`python -S research/certificate_scaling/residual_domination.py` writes
`results/certificate_scaling/residual_domination/receipt.json`.

For the asymmetric held-out edge set
`(01,1), (02,3/4), (13,1/2), (23,1/4)` at `N=2`, the naive allowance is
`5`, while the local rational row allowance is `7/2`, a `10/7` improvement.
The independent fixed-sector audit has row norm `5/2`, so the proposed
allowance passes the audit with slack.

The adverse star `(01,1), (02,1), (03,1)` has naive and proposed allowances
both equal to `6`.  The audit is `2`.  This is deliberately retained: the
rule is a sufficient local bound, not a universal improvement theorem.  A
certificate pipeline must permit such cases to fall back to another local
bound rather than claim that row grouping always helps.

The experiment is narrow.  It covers one-body residual blocks.  It does not
bound arbitrary two-body Coulomb residuals, prove that a short SOS certificate
exists, or make discovery efficient.  Extending it to density-density blocks
requires a separate local block norm bound and must preserve exact rational
verification.

## Actual molecular H4 fixture

The same regrouping was applied to every rank-one quantized block in
`results/marginal_molecule_stress/h4_square_degree3_certificate.json`.  The
full residual remains charged by exact coefficient-
`l1` for all higher-body words; only degree-two one-body words are grouped.
The particle-hole shift was included exactly: for a one-body matrix `A`, write
`A=A_0+alpha I` with `A_0` traceless and charge `alpha*N` exactly.  The
traceless part uses `min(N,M-N)` in the row bound, justified by particle-hole
duality on a traceless number-conserving one-body operator.

The total naive residual allowance is
`464099859503235409/10^18` (about `0.46410`).  The mixed bound is
`115169790838142689/2.5*10^17` (about `0.46068`), a modest exact improvement
of about `0.74%`.  Only one of 25 blocks improves; this is therefore a real
but weak result on the actual molecular fixture, not evidence that the whole
H4 residual is compressible.  The script is run with ordinary Python because
the fixture audit uses NumPy; the verifier arithmetic and receipt values are
still rational.
