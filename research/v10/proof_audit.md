# V10 proof audit: bounded independent checking of the Duhamel expansion

This note fixes the mathematical certificate that an independent checker must
replay. It assumes the v9 reference model: on the finite-dimensional Hermitian
operator space, `G = G0 + B`, `B(X) = i[V,X]` with `V=V†`, and

```text
D0(t) = exp(t G0) O,
Dk(t) = integral_0^t exp((t-s)G0) B D{k-1}(s) ds   (k >= 1),
Pm(t) = sum_{k=0}^m Dk(t).
```

The checker should validate the input assumptions and a finite, exact mode
certificate; it should then evaluate only scalar exponential bounds and the
final norm certificate. It need not build `B Dm` as a matrix merely to check
the error bound.

## 1. Differential identities and the residual certificate

Finite-dimensional differentiation under the integral gives

```text
D0' = G0 D0,       D0(0)=O,
Dk' = G0 Dk + B D{k-1},   Dk(0)=0  (k>=1).
```

Therefore

```text
Pm' - G Pm = -B Dm,       Pm(0)=O.
```

If `U(t)=exp(tG)`, variation of constants gives the exact identity

```text
U(T)O - Pm(T) = integral_0^T U(T-s) B Dm(s) ds.
```

The sign in the residual is immaterial after taking norms, but the identity is
an important audit check: a purported recurrence or coefficient convention
that produces `+B Dm` has reversed the truncation residual.

Assume the physical semigroups `U(t)` and `U0(t)` are positive, unital, and
operator-norm contractive for `t>=0`. With

```text
K >= ||B||_{infinity -> infinity},
```

induction gives

```text
||Dk(t)||infinity <= ||O||infinity (K t)^k/k!.
```

The base case uses contractivity of `U0`. For the step, use the integral
definition, contractivity of `U0(t-s)`, the induced-norm bound on `B`, and
the integral of `s^(k-1)/(k-1)!` from zero to `t`. Applying contractivity of
the *full* `U` to the residual identity yields

```text
||U(T)O-Pm(T)||infinity
 <= ||O||infinity (K T)^(m+1)/(m+1)!.
```

There is no `exp(KT)` factor. Such a factor comes from replacing the given
contractive full propagator by a generic growth estimate. Adding it is a
valid weakening only if explicitly stated; omitting it is the sharper result
under the stated hypotheses. The checker must not silently use a scalar
eigenvalue estimate for `B` or `G`: these maps may be nonnormal.

For the commutator,

```text
||B(X)||infinity = ||i[V,X]||infinity <= 2 ||V||infinity ||X||infinity,
```

so `K=2||V||infinity` is a sound, usually conservative, choice. The checker
must validate `K >= 0` and the claimed bound, rather than infer it from an
unverified spectrum.

## 2. Exact reference modes and coefficient checks

After the permitted local Clifford changes, `G0` has onsite Z fields and the
same local depolarization. In the one-qubit eigenbasis

```text
I, Z, P=(X+iY)/2, M=(X-iY)/2,
```

the local eigenvalues are respectively

```text
0,  -gamma,  -gamma+2 i h,  -gamma-2 i h.
```

Tensor products add these eigenvalues. Let the resulting exact frequencies be
`lambda_a`; then the coefficient representation is a finite sum of terms
`c exp(lambda t)` times powers of `t` generated only by repeated frequency
collisions. Frequencies and coefficients must be represented exactly (for
example, rational real and imaginary parts when that is the input domain).
Distinct terms with equal frequencies must be combined, or evaluated by the
correct confluent limit. Dividing by `lambda_a-lambda_b` when that difference
is zero is unsound.

The independent checker should verify the certificate recursively, by applying
the following exact tests to every output coefficient:

* `D0` has exactly the input expansion in the tensor-product eigenbasis,
  multiplied termwise by `exp(lambda_a t)`. At `t=0`, its coefficients equal
  the decomposition of `O`.
* For every `k>=1`, each mode coefficient in `Dk'` equals the coefficient in
  `G0 Dk + B D{k-1}`. This is an identity of exponential-polynomials, so the
  checker compares canonical frequency/polynomial coefficients exactly.
* `Dk(0)=0` for every `k>=1`, including all collision-generated polynomial
  terms. This catches integration constants and missing lower-degree terms.
* The recurrence is checked through the highest retained order, including
  `D_m`; the final norm proof still uses the bound on the complete `D_m`, not
  a guessed residual assembled from an omitted mode.

An equivalent integration check may use the exact antiderivative rule

```text
integral t^r exp(mu t) dt,
```

with a zero-at-zero constant and exact confluent handling. Recurrence checking
is preferable for an independent implementation because differentiation,
initial values, and the known action of `G0` and `B` expose errors separately.
The checker need not construct a dense `B Dm` solely to certify the final
error: it can replay the recurrence symbolically through order `m`, then use
the factorial majorant. Constructing `B Dm` is needed only if one elects to
check the residual identity numerically or symbolically as an additional
diagnostic.

## 3. Hermiticity and complete sums

For Hermitian `O`, Hermitian `V`, and adjoint-preserving `G0`, `B` preserves
Hermiticity: `G0` is a physical generator and
`(i[V,X])† = i[V,X]` when `X=X†`. Induction in the integral recurrence then
shows every complete `Dk(t)` is Hermitian for real `t`.

Individual `P` and `M` modes, and individual exponential-polynomial terms,
are complex. Hermiticity is certified only after conjugate terms are paired:
the coefficient/frequency data for a term must have its conjugate partner,
with the appropriate conjugated polynomial coefficient and the swapped
`P/M` tensor label. Self-conjugate terms must have real Hermitian operator
coefficients. Canonicalization must occur before checking this property;
approximate floating-point equality is not an adequate proof.

The factorial inequality applies to each complete Hermitian sum `Dk`, and then
to its operator norm. It cannot be applied independently to an arbitrary
unpaired complex subset of modes. A checker may instead establish a norm bound
for such a subset by a separate triangle estimate, but that is a different
certificate and must not be presented as the Hermitian contraction proof.

## 4. Input assumptions and decision boundary

The certificate is sound only after validating the parameter domain:

* `gamma >= 0`, with `gamma` real; `H` (and every local field used to form
  frequencies) is real in the declared exact domain.
* `V=V†`, `O=O†` when the Hermitian proof is invoked, and the local Clifford
  change is actually allowed by the model. A general field direction cannot
  be silently treated as a Z field by a Clifford; unsupported components must
  remain in `B` or the certificate must reject.
* `T >= 0`, finite order `m >= 0`, finite mode support, and exact dimensions
  consistent with the tensor-product basis.
* The claimed full and reference propagator contractivity assumptions are part
  of the physical model and must be preserved. The exact mode evaluator alone
  does not prove positivity, unitality, or contractivity.
* The reported `K` is nonnegative and independently justified, e.g. by
  `K=2||V||infinity` or a separately certified upper bound.

Finite support caps are a soundness boundary, not an approximation. If mode
generation, polynomial degree, coefficient size, or exact bit length exceeds a
cap, the checker must reject (or return “not certified”); it must not drop
terms, merge merely near frequencies, or substitute rounded coefficients.
Likewise, outward-rounded bounds for real exponentials and trigonometric
factors are required when numerical evaluation is used. Exact rational modes
do not make an inward-rounded transcendental evaluation rigorous.

## 5. Readout error and final composition

The Duhamel certificate bounds dynamical propagation error only. If the
experiment or protocol has an independently certified readout/measurement
error `epsilon_readout` in the same norm and at the same output interface,
the reported total error may use the triangle inequality

```text
epsilon_total <= epsilon_Duhamel + epsilon_readout.
```

The two terms must not be conflated, double-counted, or added before checking
that they refer to the same observable/state normalization. A readout bound
without an explicit conversion to the checker’s norm is not admissible.

## 6. Minimal independent replay

A bounded checker can therefore proceed as follows: validate the exact input
domain and caps; canonicalize frequencies and polynomial terms; verify `D0`,
all zero initial conditions, and the coefficient recurrence through `D_m`;
verify conjugate pairing for the complete sums when using the Hermitian proof;
validate the independent `K` and contractivity assumptions; evaluate the
factorial tail with outward rounding; and finally add a separately converted
readout bound if present. Any failed identity, unsupported parameter, omitted
collision term, cap exceedance, or unverifiable physical assumption is a
rejection, not an invitation to fall back to an unchecked approximation.

## 7. Root review of implementation audit

The agent's initial dense comparison incorrectly compared `Reference.column`,
which implements `i[V,.]`, with `i[H,.]` including the removed onsite field.
That mismatch does not establish a conversion defect. `expand_word` and the
local-product expansion concatenate one character per tensor factor, so each
output path is unique within that operation. Contributions from different
input terms are accumulated explicitly.

Root corrected the comparison to V and enumerated all 4^n eigenwords for n1
and n2, rather than only the words present in V. Those direct dense comparisons
pass. An exact three-qubit Pauli/eigenbasis round trip also passes. Two initial
agent tests had invalid construction expectations (an unattainable tolerance
at order3, and a request for a D1 layer when V=0). They were replaced by real
recurrence, initial-condition and Hermiticity mutation tests. The original
failed audit/tests are preserved in `invalidated/`; they are not evidence of
a defect in the surviving algebra or proof.
