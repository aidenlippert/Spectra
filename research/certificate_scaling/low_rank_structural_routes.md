# Structural routes beyond raw Gram rank

This note separates four compression ideas that are often conflated:

* Coulomb or density fitting controls the number of interaction integrals. It
  does not control the rank of the positive SOS Gram after CAR normal ordering.
* Tensor/MPO or operator-Schmidt rank controls a representation across a cut.
  It does not imply a small global PSD factor, a small spectral gap, or an
  easy ground-state search.
* Basis adaptation can make a particular state or operator compressible, but
  the transformation must be included in the certificate and its truncation
  error bounded in the original operator norm.
* Exact algebraic squares are the strongest route: if the Hamiltonian is
  explicitly a sum of a small number of squares modulo an exact ideal, the
  certificate is already discovered and the proof is cheap.

## A checkable sufficient condition

Let (H) be represented in a fixed CAR word basis and suppose an exact
decomposition is

\[
 H-b = \sum_{a=1}^r Q_a^\dagger Q_a + (\hat N-N)X + R.
\]

If the retained factors are rational and the omitted operator satisfies a
checkable induced norm bound \(\|R\|\leq\varepsilon\), then

\[
 E_0(H)\geq b-\varepsilon.
\]

For a normal-ordered expansion (R=\sum_w c_w w), the conservative bound
\(\varepsilon=\sum_w|c_w|\) is immediately checkable because every CAR
monomial has operator norm at most one. A weighted basis can improve this to
\(\sum_w |c_w|\,\|w\|_{\rm bound}\), but the weights must themselves be
proved. Thus a useful structural theorem must give, uniformly in system size,

1. a discoverable factor count (r(M)),
2. an exact or outward-rounded residual representation, and
3. a tail bound with 
   \(\sum_w|c_w|\leq\varepsilon(M)\) at the desired total energy accuracy.

Tensor rank or integral sparsity alone supplies none of these three items.

## Exact overlapping quartic family

For an even periodic chain of (M) spinless modes, define

\[
 q_i=n_i+n_{i+1}-1,\qquad H_M=\sum_{i=0}^{M-1}q_i^2,
 \quad n_i=a_i^\dagger a_i,
\]

with indices modulo (M). Each (q_i) is a number-conserving quadratic
operator, so (q_i^2) contains genuine quartic density terms and overlapping
neighbours share modes. The factors have rank one per bond and total storage
is (O(M)), while the interaction graph is connected rather than a union of
independent blocks.

The alternating Fock state
\(|1010\ldots\rangle\) has (q_i|\psi\rangle=0) on every bond. Therefore
the exact SOS lower bound (b=0) and this explicit upper witness both give
\(E_0=0\) in the half-filled sector. Discovery is scalable: emit the same
local factor for each translated bond, and the exact replay checks each local
CAR identity plus the periodic closure. No Hilbert-space enumeration is
needed.

This family demonstrates the strongest positive case for structural
compression, but it is intentionally frustration-free. Adding a hopping term
or changing one bond generally destroys the common kernel; the same local
factors still prove a lower bound only after paying the exact residual norm.
The family therefore does not imply that generic molecular or Hubbard
Hamiltonians have (O(M)) certificates, nor that a low MPO rank solves their
ground-state search.

## Matched upper bounds and total accuracy

For a nonzero tail, a matched upper bound must be computed from an explicit
state or variational ansatz and independently certified. The meaningful total
interval is

\[
 U_M-(b-\varepsilon)= (U_M-b)+\varepsilon.
\]

Controlling only the SOS tail leaves (U_M-b) uncontrolled. Conversely, a
small variational error does not repair an uncertified lower bound. A route is
therefore promising only when both terms have a uniform scaling law and the
factor discovery cost is included.

The practical next test is to perturb the overlapping family by bounded
finite-range hopping, generate exact rational residuals for retained local
squares, and measure whether the induced (L^1) tail grows linearly, remains
bounded, or destroys the useful interval. That experiment distinguishes a
real structural theorem from a frustration-free special case.

## Executable noncommuting control

`research/certificate_scaling/low_rank_structural_check.py` implements a
genuinely overlapping family using two modes per site,

\[
B_i=a_{\downarrow,i}a_{\uparrow,i+1}
    -a_{\uparrow,i}a_{\downarrow,i+1},\qquad H=\sum_i B_i^\dagger B_i.
\]

Adjacent factors share both modes on their common site, so the local factors
are not disjoint. The all-up determinant at particle number (N=L) is
annihilated by every (B_i), giving an exact lower and matched upper bound of
zero. The standard-library generator and repository exact verifier run for
(L=4,8,16,32), with factor storage (2L) and no source factors or Fock
matrix. This is a structured control family, not a claim about generic
chemistry; its usefulness is that it tests whether a proposed compression
argument survives overlapping noncommuting factors before being applied to
hard Hamiltonians.
