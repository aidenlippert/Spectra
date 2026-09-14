# Dual signs and a safe pricing rule

## LP and dual

Write the coefficient-space relaxation as

```
maximize    b - 1^T t
subject to  h = b e_0 + C lambda + A x + r
            -t <= r <= t,   lambda >= 0, t >= 0,
```

with `b,r,x` free. Here columns of `C` are normal-ordered SOS atoms, columns of
`A` are explicitly supplied number/CAR ideal maps, and `h` is
the Hamiltonian coefficient vector. The primal objective is the lower endpoint
`b-||r||_1`; it is only a physical lower bound when the retained identity and
residual have a separate positivity/operator-norm verification.

Using a multiplier `y` for `h-b e_0-C lambda-r`, the dual is

```
minimize    h^T y
subject to  e_0^T y = 1,
            C^T y >= 0,
            A^T y = 0,
            ||y||_infinity <= 1.
```

The sign is easy to test: a new atom column `c` has dual slack `c^T y`. Under
the repository convention `minLPdual=-y`, the same reduced-cost test is written
`-c^T minLPdual < 0`; convert once at the interface and do not mix conventions.

## Operator-word pricing

For words `w_i`, an atom `B=sum_i z_i w_i` has coefficient column `c(z)` from
the exact CAR normal-form map of `B^dagger B`. The reduced cost is

```
c(z)^T y = z^dagger Q(y) z,
Q_ij(y) = y^T coeff(w_i^dagger w_j).
```

Therefore a negative eigenvalue of the Hermitian `Q(y)` supplies a violating
column: its normalized eigenvector `z` gives `z^dagger Q z<0`. In exact work,
an interval or rational Rayleigh quotient must certify strict negativity after
CAR expansion; a floating negative eigenvalue is only a pricing hint.

This guarantees a dual-feasibility violation, not a guaranteed energy gain of a
specified size. The restricted LP must be re-solved, and the resulting primal
certificate must pass its independent positivity and residual checks.

## Omitted-family bound with a trace budget

Suppose all omitted atoms are indexed by `a`, with PSD weights `lambda_a`, and
the search imposes a known factor budget

```
sum_a lambda_a <= B.
```

If a dual point has certified quadratic forms

```
z^dagger Q_a(y) z >= -delta
```

for every normalized omitted atom, then every omitted SOS contribution has dual
violation at most `delta`, and the total omitted contribution is bounded by
`delta B`. Thus an omitted-family pricing certificate can safely report an
additional objective uncertainty `delta B`, provided the same normalization and
budget are enforced by the primal. A per-atom bound without a global budget is
not enough when the family is unbounded.

### Proof

The dual pairing of an omitted weighted atom is
`lambda_a z_a^dagger Q_a z_a >= -lambda_a delta`. Summing and using the budget
gives a total of at least `-delta sum_a lambda_a >= -delta B`. This is the only
step that converts a uniform per-normalized-atom bound into a family bound.

## Why a small negative eigenvalue is not a small objective certificate

Take a two-word pricing matrix

```
Q = [[-epsilon, 0], [0, 1]],   epsilon > 0,
```

whose normalized negative eigenvector is `(1,0)` and reduced cost is
`-epsilon`. This by itself is only a pairing statement. The fully explicit
**abstract coefficient LP**

```
maximize b-|r|,       0 = b - epsilon*lambda + r,   lambda >= 0
```

has dual test value `-epsilon` for the new column and is unbounded above:
choose `r=0`, `b=epsilon*lambda`, and let `lambda` grow. Thus a small negative
reduced cost does not imply a small objective effect when no weight bound
exists. If the same LP adds `0 <= lambda <= B`, its optimum is exactly
`epsilon B`. This is not a physically valid SOS atom by itself: it is a
coefficient-space counterexample showing why reduced cost alone cannot
quantify objective gain. Scaling the atom changes both its coefficient column and the
meaning of its weight, so normalization must be fixed before comparing
eigenvalues. With an explicit trace budget the total improvement is bounded by
`B epsilon` under that fixed normalization.

This is the exact reason pricing must record word normalization and a factor or
trace budget. A negative eigenvalue ranks a direction; it does not quantify the
final bound improvement, omitted-family error, or physical certificate quality.

## Practical interface

1. Export the LP dual in one declared sign convention.
2. Build `Q(y)` by exact CAR coefficient maps, preserving the word normalization.
3. If the direction is to be claimed as a *certified pricing violation*, require
   interval/rational certification that its Rayleigh quotient is negative. A
   heuristic direction may still be passed to the primal solver without that
   proof, but it earns no pricing guarantee; only a subsequent exact SOS, ideal,
   and residual replay can make the resulting certificate valid.
4. Re-solve the restricted LP and independently verify the SOS identity,
   sector-safe number relation, and residual operator bound.
5. Use `delta B` only when a finite factor budget and a normalized omitted-family
   bound are both explicit.

No step above proves that the dual-guided dictionary search is polynomial or
that generic molecular Hamiltonians have a compact certificate.
