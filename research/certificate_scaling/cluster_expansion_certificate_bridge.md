# Conditional bridge from cluster tails to fixed-total-error certificates

This note isolates the extra statement needed to turn a cluster expansion into
a positive SOS lower certificate.  Local energy approximation alone is not
enough: a numerical or perturbative approximation to `E0` does not provide a
one-sided operator inequality.

## Conditional tail theorem

Let a finite system have `N` sites and a decomposition

`H = H_(<=k) + R_(>k)`.

Assume every connected cluster `X` of size `ell` in the tail has a certified
operator decomposition

`R_(>k) = sum_X R_X`,   `||R_X|| <= A exp(-mu ell)`,

and assume the number of connected clusters of size `ell` is at most
`N C^(ell-1)`.  Then, provided `mu > log C`,

`||R_(>k)|| <= N A exp(-(mu-log C)(k+1)) / (1-exp(-(mu-log C)))`.

This is just the triangle inequality plus the cluster count, but it is a
useful exact accounting identity.  If the truncated Hamiltonian has an exact
positive certificate

`H_(<=k) - b_k I = S_k + J_k`,   `S_k >= 0`,

then the full Hamiltonian has the certified lower bound

`E0(H) >= b_k - eta_k`,

with `eta_k` equal to the displayed tail bound.  The verifier need only check
the local rational inequalities for each `R_X`, the exact SOS replay for
`H_(<=k)`, and the geometric-series arithmetic.

To obtain total additive error at most `epsilon` from the tail, it is enough
to choose

`k >= log(A N / (epsilon (1-exp(-nu)))) / nu - 1`,

where `nu = mu-log C`.  Thus `k=O(log(N/epsilon))` for fixed `A,mu,C`.

## Conditional size and bit-cost accounting

Suppose the number of retained connected terms of size at most `k` is at most
`N C^k`, each local term uses at most `s(k)` rational coefficients, and every
coefficient has `B(k,p)` bits when the target accuracy is `p` bits.  Direct
local expansion and rational checking then cost

`O(N C^k s(k) poly(k) B(k,p))`

arithmetic operations and the certificate has the same term count up to the
SOS factors.  Substituting the logarithmic `k` gives a polynomial in
`N/epsilon` for fixed `nu`, `C`, and polynomial `s`.  The exponent is explicit:
`C^k = O((A N/epsilon)^(log C/nu))`.  Rational geometric tails require only
`O(log(N/epsilon))` additional bits when `A,C,mu` are supplied as certified
rational enclosures.  These bounds count coefficient discovery separately;
they do not make finding the local factors efficient.

## The missing bridge to a positive SOS certificate

A ground-state cluster expansion does not imply the displayed identity.  A
sufficient stronger hypothesis is a certified quasilocal similarity/unitary
normal form: there are local generators `G_X` and a truncated transform
`U_k=exp(G_(<=k))` such that

`U_k H U_k^dagger = E_k I + sum_a Q_a^dagger Q_a + R_(>k) + J_N`,

where every `Q_a` and every omitted cluster has a local rational description,
`J_N` is the fixed-number ideal, and the tail obeys the certified bound above.
Conjugating a positive square preserves positivity; the only one-sided error is
the explicitly bounded tail.  A nonunitary similarity transform needs an
additional condition-number bound and cannot be treated as positivity
preserving for free.

The existence of such a quasilocal positive normal form is **not proved here**.
Neither a ground-state cluster expansion nor a generic correlation-length
claim supplies it.  In particular, no generic gapped assumption is being
smuggled in: a gap, locality, screening, and bounds on nested commutators may
be useful hypotheses for proving the normal form, but each must be stated and
certified.  Gapless systems, critical points, long-range Coulomb couplings,
and basis changes can violate the assumed exponential tail.

The falsifiable research target is therefore precise: construct the local
positive normal form and prove its rational tail bounds for a nontrivial
chemical family.  Until that construction exists, the theorem above is a
conditional scaling result, not a claim that cluster expansions already solve
certificate discovery.

## Verification boundary

The uniform omitted-cluster estimate must itself follow from efficiently
checkable Hamiltonian hypotheses or an analytic tail certificate. Enumerating
and checking every omitted cluster would defeat the stated complexity bound.
Likewise, conjugating squares by a quasilocal unitary preserves positivity,
but does not automatically preserve compact factor storage; truncation,
coefficient precision and expansion costs for the pulled-back factors require
their own bounds. A nonunitary similarity transformation is not a congruence
and does not preserve a Hermitian PSD ordering merely because its condition
number is bounded. A valid metric/congruence argument would be extra work.
