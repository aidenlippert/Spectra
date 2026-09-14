# Density patterns with a collectively bounded remainder

All claims concern the supplied rational electronic Hamiltonian on its fixed
total-particle-number sector. The molecular fixtures use interleaved spin
orbitals. No spin-purity assumption or many-body eigenvector is needed for
acceptance. The implementation is `core.py`.

## 1. Reconstruct the interaction before compressing it

For s spatial orbitals let

\[
E_{pq}=\sum_{\sigma=\uparrow,\downarrow}a^\dagger_{p\sigma}a_{q\sigma}.
\]

Use the d=s(s+1)/2 real Hermitian features Q_pp=E_pp and
Q_pq=E_pq+E_qp for p<q. Opposite-spin quartic coefficients give an electron
repulsion matrix C on these features; the extraction averages the internal
pair permutations and pair exchange. Rather than assuming a convention for
one-body contraction terms, expand ½Σ_ab C_ab Q_a Q_b using exact CAR.
Subtract its constant and quadratic coefficients from the original H to
define h_one. Reconstruct the remaining quartic mismatch W exactly. Thus

\[
H=h_{\rm one}+\tfrac12 Q^T C Q+W,\qquad \|W\|\le\eta. \tag{1}
\]

The norm bound sums absolute diagonal-monomial coefficients and one magnitude
per conjugate off-diagonal pair. A diagonal CAR monomial is a signed occupation
projector. An off-diagonal number-conserving monomial maps disjoint occupation
subspaces; the norm of a Hermitian conjugate pair is at most its coefficient
magnitude. This supplies a full-Fock-space sufficient bound.

For the frozen H6 data η=3.2×10⁻¹¹ Ha. It measures reconstruction mismatch to
the rational input, not basis-set error or all error in an experimental model.
C is not assumed exactly PSD: tiny negative values from rational coefficient
rounding are covered by the subsequent exact envelope.

## 2. Eliminate the conserved number direction exactly

Let u_a=1 on diagonal features and zero elsewhere, so uᵀQ=N̂ and uᵀu=s.
Set P=I−uuᵀ/s and Q₀=PQ. On the N̂=n sector,

\[
\tfrac12 Q^TCQ
=\tfrac12 Q_0^TCQ_0
+\frac n s (PCu)^T Q
+\frac{n^2}{2s^2}u^TCu. \tag{2}
\]

The last two terms are respectively one-body and constant terms. Fold them
into h̃_one and use C₀=PCP for factor discovery. This is an operator identity
on the fixed-number sector, including states with indefinite individual
orbital occupations. It preserves every coupling to the number direction.

The program also runs an uncentered control using (1) directly. A centered
factor vector v is required to satisfy uᵀv=0 **exactly**, including after
rational rounding. Otherwise the improved bound below would not apply.

## 3. Collectively bound the remaining features

Let G be diagonal with G_aa=1 for diagonal features and 1/2 for off-diagonal
features. The orbital and spin Casimir identity is

\[
\sum_{pq}E_{pq}E_{qp}
=(s+2)\hat N-\tfrac12\hat N^2-2\mathbf S^2. \tag{3}
\]

One derivation expands the orbital expression and the two-spin-species
expression Σ_στ X_στ X_τσ, where X_στ=Σ_p a†_pσ a_pτ. Their quartic terms
cancel, giving (s+2)N̂; the second expression equals N̂²/2+2S².

For p<q let T_pq=E_pq−E_qp. Separating real symmetric and antisymmetric
features gives

\[
\sum_{pq}E_{pq}E_{qp}
=Q^TGQ+\tfrac12\sum_{p<q}T_{pq}^\dagger T_{pq}.
\]

Also Q₀ᵀGQ₀=QᵀGQ−N̂²/s. Positivity of S² and every T†T therefore proves,
on the fixed-n sector,

\[
0\preceq Q_0^TGQ_0\preceq
\kappa I,\qquad
\kappa=n(s+2)-n^2/2-n^2/s. \tag{4}
\]

For uncentered features omit the last n²/s subtraction. The identities were
also checked by independent exact CAR expansions on two and three spatial
orbitals. Equation (4) bounds all complementary occupations collectively;
its proof does not enumerate the physical sector.

## 4. Rational factor and remainder acceptance

The retained coefficient matrix is K=Σ_{k=1}^r λ_k v_kv_kᵀ with rational
λ_k≥0. Each retained operator A_k=v_kᵀQ is a spatially extended density
pattern. H_r=h̃_one+½Σ_k λ_k A_k². For R=C₀−K, accept rational α≤0≤β only if

\[
R-\alpha G\succeq0,\qquad \beta G-R\succeq0. \tag{5}
\]

PSD of a coefficient matrix implies PSD of the corresponding operator form
even though the Q_a do not commute: factor the matrix as BᵀB, giving a sum
of operators (Σ_a B_ka Q_a)†(Σ_b B_kb Q_b). This argument keeps the full
operator products. It is not a one-particle positivity assertion about an
unrelated interacting occupation Hamiltonian.

Since R=PRP, combine (1), (4), and (5) to obtain

\[
H_r+\ell I\preceq H\preceq H_r+uI,
\quad
\ell=\alpha\kappa/2-\eta,\quad
u=\beta\kappa/2+\eta. \tag{6}
\]

The negative sign of α matters: α Q₀ᵀGQ₀ ≥ ακI. The certified remainder
width is (β−α)κ/2+2η. Equation (6) bounds every state and hence every ordered
eigenvalue in the sector. It transfers any valid retained-model energy
interval [L_r,U_r] to [L_r+ℓ,U_r+u]. It does not itself solve H_r.

Numerical diagonalization of G⁻¹ᐟ² C₀ G⁻¹ᐟ² proposes eigenfactors. Vectors
and nonnegative weights are rounded rationally, exact zero trace is restored,
and α/β are adjusted until (5) passes rational symmetric elimination. The
checker rejects negative pivots and nonzero coupling at a zero pivot. It
reconstructs the original input, contractions, mismatch, number identity,
factor matrix, metric, and both PSD checks; no saved floating eigenvalue can
accept a certificate.

## 5. Exact minimum-rank obstruction for this remainder bound

For a desired width ε>2η define τ=2(ε−2η)/κ. A rational d×r matrix Z with

\[
Z^T(C_0-\tau G)Z\succ0 \tag{7}
\]

proves that any retained K of rank at most r−1 fails this width budget in
the envelope (5)–(6). Indeed, its kernel intersects the r-dimensional image
of Z in a nonzero vector z. Then zᵀKz=0 and zᵀC₀z>τ zᵀGz, forcing β>τ
whenever βG−(C₀−K) is PSD. With α≤0, the certified width exceeds ε.

Exact PSD elimination with rank r verifies strict positivity in (7).
Combining this obstruction with a passing r-factor construction establishes
the minimum rank **for this coefficient envelope and collective scalar
bound**, including arbitrary orientations of the retained rank-r matrix.
It does not exclude better operator inequalities, other collective identities,
different basis sets, or different pattern languages.

## 6. A cheap retained-model solver and its exact limitation

Write h̃_one=c_one+dΓ(t). For any real centers c_k, positivity of
(A_k−c_k)² gives

\[
H_r\succeq \widetilde h_{\rm one}+\sum_k\lambda_k c_k A_k
-\tfrac12\sum_k\lambda_k c_k^2. \tag{8}
\]

The right side is one-body plus a constant. Let its one-body matrix be B.
A certificate supplies μ and X=VVᵀ+δI with δ≥0, and checks B−μI+X PSD.
Fermionic occupations lie between zero and one, so dΓ(X)≤tr(X)I. This
proves the lower endpoint

\[
b=c_{\rm one}-\tfrac12\sum_k\lambda_kc_k^2
+n\mu-\operatorname{tr}X+\ell. \tag{9}
\]

The exact upper witness used in these new intervals is the HF determinant.
Previously computed many-body reference uppers are kept separate.

To test whether weak lower bounds are a search failure, supply a rational
one-body matrix γ satisfying 0≤γ≤I and tr γ=n. For every B, the sum of its
n smallest eigenvalues is at most tr(Bγ). Maximizing the resulting scalar
quadratic expression over all centers in (8) gives the exact ceiling

\[
\sup_{c}\ b_{\rm ideal}(c)
\le c_{\rm one}+\operatorname{tr}(t\gamma)
+\tfrac12\sum_k\lambda_k[\operatorname{tr}(A_k\gamma)]^2+\ell. \tag{10}
\]

Here matrix notation denotes the one-body part; its constant is written
separately. The checker verifies both PSD constraints on γ and its exact
trace. The numerical convex program only proposes γ and the centers.

Finally, an independently replayed existing lower certificate L_* for the
same H proves E₀(H)≥L_*. When L_* exceeds the ceiling in (10), their exact
difference is a lower bound on the gap between the true ground energy and
**every** lower bound in this supporting-plane family. This is the reported
dual obstruction. It leaves joint fermionic constraints between patterns
fully open.

## 7. Resource and physical scope

The coefficient PSD checks have dimension s(s+1)/2. Rank witnesses have
dimension r, and supporting-bound PSD checks have dimension 2s. A factor has
up to s(s+1)/2 rational coefficients; factor count alone does not describe
its cost. Reading the original Hamiltonian, CAR reconstruction, exact number
folding, factor discovery, and replay are all charged. The implementation
refuses more than ten spatial orbitals in this pass.

No accepting or discovery path builds a physical many-body basis. Independent
full-sector controls are capped at twelve spin orbitals. The numerical
ground-energy changes in those controls are distinct from (6). H4–H10 form
four finite-size observations at one geometry/basis prescription, not an
asymptotic theorem or an accuracy guarantee for experimental molecules.
