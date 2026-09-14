# Composing the response without a determinant matrix

## Claim and scope

At the preserved frozen-H6 target b = U − 0.001 Ha, a second eliminated
sector of the first retained operator has a strictly positive, exactly
replayed lower bound. Its response and uniform residual bound are constructed
from rational orbital matrices and canonical fermionic operator words.
Neither construction nor acceptance enumerates the complete particle sector.
The terminal retained operator is **not** yet certified positive.

This uses established Schur-complement/Feshbach elimination machinery.
The experiment concerns its constructible operator bounds and composition,
not novelty of the elimination identity. See
[Dusson, Sigal and Stamm (2021)](https://arxiv.org/abs/2105.02058).

## First elimination and the exact remaining obligation

Let Q1 project onto double occupation of the last spatial orbital, P1 = I−Q1,
and write H−b in P1,Q1 order as [[A,B],[B*,D]]. The inherited certificate
proves δ1 I ≤ D ≤ M1 I. Define

    r_k(d) = T_k((c−d)/a) / T_k(c/a),
    p_k(d) = (1−r_k(d))/d,  c=(M1+δ1)/2, a=(M1−δ1)/2,
    X1 = p_k(D) B*,
    K1 = A − B X1 − X1* B* + X1* D X1
       = A − B F1(D) B*,  F1(d)=2p_k(d)−d p_k(d)^2.

The quotient defining p is a polynomial of degree k−1. For d in the certified
interval, |r_k(d)| ≤ 1/T_k(c/a) < 1, hence

    0 ≤ F1(D) = (I−r_k(D)^2) D^−1 ≤ D^−1.

With residual R1=B*−DX1, the exact Schur complement is

    S1 = K1 − R1* D^−1 R1 ≥ K1 − η1 I = L1,
    η1 = g1^2 / (δ1 T_k(c/a)^2).

Consequently L1≥0 is sufficient for H≥b. No eigenvector or matrix of L1 is
needed to state or verify the uniform residual estimate.

## Orbital-sized filling certificate

The input tail certificate supplies an exact operator sandwich around

    c0 + dΓ(t) + (1/2) Σ_l w_l dΓ(L_l)^2,  w_l>0.

The unretained density remainder is handled collectively by its certified
lower/upper shifts. Each spatial matrix acts identically on two spin species.
For a rational spatial matrix S, rational τ and rational symmetric Y with

    Y≥0,   S−τI+Y≥0,

fermionic occupancy at most two per spatial orbital gives

    dΓ(S) ≥ τ N − dΓ(Y) ≥ τ N − 2 tr(Y).

The accepting checker verifies the two small PSD matrices with exact
rational arithmetic. Numerical optimization only proposes tangent values,
τ and Y; no numerical eigenvalue is trusted as a lower bound.

For a single doubly occupied orbital j, write L=[[A,v],[v*,ell]]. Compression
to that occupation sector gives the exact identity

    Qj dΓ(L)^2 Qj
      = Qj[(dΓ(A)+2ell)^2 + 2||v||² − dΓ(vv*)]Qj.

For each rational tangent x, use
(dΓ(A)+2ell)^2 ≥ 2x(dΓ(A)+2ell)−x², then apply the filling certificate.
This strengthens the first denominator bound without a determinant basis.

## Why a separate coupling norm loses too much information

Let Q2 be double occupation of the next spatial orbital and R2=P1 Q2.
Writing α2 for a bare bound on Q2(H−b)Q2, a simple sufficient bound is

    R2 L1 R2 ≥ (α2 − ν/δ1 − η1) R2,
    ν ≥ ||Q1 H R2||².

`coupling.py` obtains ν using seven occupation labels on only four local spin
modes. The remaining modes stay as CAR operators. Odd transition blocks
use {C*,C} to cancel the highest degree terms; even blocks use C*C. Exact
coefficient bounds, outward square roots and block row/column norm bounds
complete the proof. H6's largest local majorant has 184 words, H8's 919.

This estimate is negative on both fixtures. That alone does not prove the
actual second-sector gap absent. `recursive.py` additionally constructs
an exact two-dimensional example [[a,t],[t,d]] consistent with the supplied
α2, δ1, M1 and coupling-norm bounds but with

    a − t²(2p_k(d)−d p_k(d)^2) < 0.

Thus those separate scalar facts cannot force the required positivity.
This is an obstruction to sufficiency of the scalar information, not a
family optimality witness or a counterexample to the molecular method.

## Joint occupation bound: carry the induced term intact

Set R = Q1 ∨ Q2 = Q1+R2. Every vector in R has at least two electrons in
the last two spatial orbitals. Let P_high be their one-particle projector.
For rational x_l and μ≥0, density-square tangents give the global inequality

    H ≥ c0 − (1/2)Σ_l w_l x_l²
        + dΓ(t+Σ_l w_l x_l L_l−μ P_high)
        + μ N_high + tail_lower.

On R, replace μ N_high by 2μ and apply the small filling certificate. This
proves R H R ≥ γ_joint R, using a 6×6 orbital matrix for H6 and 8×8 for H8.
The shared tail replay additionally checks coefficient matrices of dimensions
21 and 36 respectively; these also scale with orbitals, not determinants.
The constraint N_high≥2 is a relaxation of the union-of-double-occupation
sector, so the bound is valid but need not be sharp.

**Joint-sector lemma.** If R(H−b)R ≥ γ R with γ>0, then for any response
X1 mapping P1 to Q1,

    R2 K1 R2 ≥ γ R2.

Indeed, for v in R2, the vector (v,−D^−1 B*v) belongs to R and has norm at
least ||v||. Its H−b quadratic form is v*S1v, so the restricted exact
Schur complement is at least γ. The identity

    K1 = S1 + (X1−D^−1B*)* D (X1−D^−1B*)

then proves the claim. The inverse is used in the argument; the implemented
certificate constructs no inverse or determinant matrix.

It follows directly that

    D2 := R2 L1 R2 ≥ (γ_joint−b−η1) R2.

For frozen H6 the exactly checked value is **0.04023822566844535 Ha**.
This avoids replacing the induced interaction by ν/δ1 in the gap proof.

## Second response and its remaining cost

Let P2=P1(I−Q2). Since F1≥0, D2 is at most the corresponding compression of
H−b−η1; the existing density-square upper estimate provides M2. For the
off-diagonal coupling of L1 one may still use

    g2 ≤ ||Q2 H(I−Q2)||
         + ||Q1 H R2|| ||Q1 H P2|| / δ1_sharp.

Only the response error estimate uses this norm bound. δ1_sharp comes from
the exact single-orbital tangent certificate. The first response still uses
its original certified spectral interval; 0≤F1≤D^−1 remains valid.

Outward-rounded H6 parameters are δ2=0.0402 Ha, M2=14.4518 Ha and
g2=28.2965 Ha. The same Chebyshev stopping rule gives k2=119 and
η2≈9.7328641091e−7 Ha. Apply the preceding Schur argument to L1, using
X2=p_119(D2) times its coupling. This leaves the sufficient obligation
K2−η2 P2≥0. Its proof has **not** been constructed.

`recursive.retained_action` expresses K through parent actions and occupation
projections. The exact two-level action/congruence identity is tested on a
small rational example. No molecular second response was expanded into
vectors or evaluated on all columns. The implemented recurrence takes
2k1+1=53 parent H actions per K1 action. Therefore an outer polynomial
evaluation on a *supplied* right-hand side costs 118×53=6,254 H actions;
the complete K2 action costs 239×53=12,667 H actions. Obtaining the coupling
right-hand side is additional to the former count and included in the latter.
These are counts for this straightforward oracle implementation, not a
lower bound on every possible implementation. Intermediate vector support,
rational bit growth and actual action cost remain unresolved.

Combinatorially, P1 has dimension 714 and P2 dimension
924−210−(210−28)=532 on frozen H6. These counts are not constructed bases.
The result is a composable, enumeration-free second response component;
the full retained proof still has an open terminal obligation.
