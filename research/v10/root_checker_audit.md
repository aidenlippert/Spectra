# v10 readout/root checker audit

Scope: read-only inspection of `experiments/v10_reference.py` and
`experiments/v10_exp_readout.py` as currently written. This is an audit of
proof conditions and failure modes, not a global certification from tests.

## Findings

### 1. The exponential disk arithmetic is sound under its stated l1 norm

`norm(x+iy)=|x|+|y|` is submultiplicative, so the disk product radius
`||a|| r + ||b|| s + rs` is valid. In `exp_disk`, after scaling
`w=z/2**scalings`, `||w|| <= 1/2`; for Taylor degree `n`, the next-term ratio
is `||w||/(n+1) < 1`. Therefore

```
||term_n|| * ratio/(1-ratio)
```

is a valid geometric upper bound for the omitted tail. (The code starts at
`n=1`, and the ratio is at most 1/4.) The midpoint rounding adds its exact l1
distance to the radius. Each squaring uses
`2*||mid||*radius + radius**2`, which is the correct enclosure for squaring a
complex disk. The repeated grid rounding is also charged. The scaling loop is
bounded and refuses unsupported requests.

The use of floor for each signed grid coordinate is conservative because the
actual midpoint-to-grid distance is added to the radius; it is not a bug.

### 2. Aggregation and exponent-weight allocation are conservative

Rows are aggregated by exact `(a,b)` and then by word, with the polynomial
factor `t**k` included before the coefficient norm is computed. Exact
cancellation removes zero coefficients. For each nonempty exponent group,
`weight = sum_word ||coefficient||_1`; the target is
`tolerance / max(1,N) / max(1,weight)`, where `N` is the number of nonempty
exponent groups. Hence the charged group error is at most
`weight*target <= tolerance/N`, and summing groups is at most the requested
tolerance. The final internal check catches arithmetic or allocation failure.

No collision-overwrite issue was found in `expand_word`: its generated keys
are concatenated local paths, and distinct paths can collide only when the
same key is intentionally combined. In the readout aggregation, `setdefault`
and explicit addition likewise preserve collisions.

### 3. The recurrence checker and proposer integration agree

`integrate` solves the coefficient ODE
`Q' + (lambda-mu)Q = c t^k`, including the `lambda == mu` degree-raising
case. `derivative_minus_reference` computes exactly the same identity and
`check_expansion` compares it against `apply_b(previous_layer)`. Exact rational
spot checks with nontrivial Hamiltonian interaction and repeated frequencies
returned equality through several recurrence levels. This supports the local
integration identity; it does not certify arbitrary caller inputs beyond the
checker’s explicit validations.

### 4. Basis conversion and real-part projection are justified at operator norm

The unchanged `error` returned by `evaluate_expansion` is sufficient for the
operator-norm interface. If the eigenbasis midpoint error is
`E = sum_W delta_c_W W`, every tensor eigenword `W` in `I/Z/P/M` has operator
norm 1. Therefore `||E||_infty <= sum_W |delta_c_W|`, and the latter is exactly
bounded by the complex l1 error accumulated by `evaluate_mode_map`. The exact
eigenoperator-to-Pauli conversion is an identity of matrices, so its internal
path count is cost bookkeeping rather than norm amplification.

For the checked Hermitian exact result `M`, converting to Pauli coefficients
and retaining `c[0]` produces `(A+A†)/2`: Pauli strings are Hermitian and the
coefficient real part is the Hermitian projection. Consequently
`||(A+A†)/2 - M||_infty <= ||A-M||_infty`, since `M=M†`. No additional
projection charge is needed. This resolves the prior suspected interface gap;
the argument depends on the output claim being operator norm at the same
matrix interface and on the checker’s Hermiticity premise.

### 5. The readout contract depends on the preceding Hermiticity check

`check_expansion` verifies conjugate pairing with swapped P/M labels and
negated imaginary exponent in every layer. That is enough to establish the
formal Hermiticity of the complete exact layer under the stated basis
conventions. `evaluate_expansion` then applies the Hermitian projection by
taking real Pauli coefficients, as proved above. Its internal scalar evaluator
does not need to enforce conjugate pairing on rounded coefficients. The
caller must establish the exact Hermitian target first; `v10_headroom.attempt`
and `replay_reference` both call the expansion checker before this readout.

## Conditional conclusion

The Taylor, scaling/squaring, rounding, aggregation, recurrence, basis
conversion, and Hermitian projection formulas are locally conservative under
the stated complex l1 and operator-norm argument. The readout is usable when
the root certificate explicitly claims operator-norm error at the returned
matrix/Pauli interface and relies on the checked Hermiticity premise. The
exact tests above are spot checks and must not be treated as held-out or
global certification.
