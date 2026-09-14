# Structural attacks on the v7 residual certificate

This note proposes two bounded certificate mechanisms that exploit structure already
present in the exact residual calculation. They are candidate evaluator arms, not
claims that a learner has found a method. Every inequality below is independently
checkable from rational Pauli coefficients; if a side condition is not proved, the
checker must fall back to the current unweighted certificate or refuse.

## 1. Bernstein cancellation for a polynomial residual

On a piece of width `h`, write a scalar residual coordinate as

```text
q(u) = sum_{k=0}^d a_k u^k,       0 <= u <= h.
```

Put `x=u/h` and convert exactly to the degree-`d` Bernstein basis,

```text
q(hx) = sum_{j=0}^d b_j B_{j,d}(x),
B_{j,d}(x) = binom(d,j)x^j(1-x)^(d-j),
b_j = sum_{k=0}^j a_k h^k binom(j,k)/binom(d,k).
```

The `B`'s are nonnegative and sum to one. Therefore the checker has the rational
witnesses

```text
sup_[0,h] |q| <= max_j |b_j|,
integral_[0,h] |q| du <= h/(d+1) sum_j |b_j|.                 (B)
```

The second bound follows by integrating the triangle envelope and
`integral B_(j,d) dx = 1/(d+1)`. For a vector of coefficients in one verified
anticommuting Pauli group, use `b_j` as vectors and replace `|b_j|` by an exact
outward square-root upper bound on `sum_i b_(i,j)^2`; the same formula applies.
For a leftover Pauli coordinate, (B) is directly applicable. Subdivision at a
rational point repeats the conversion on each subinterval and can only tighten
the certificate.

The witness must include `d,h,a,b`, with the checker recomputing `b`; supplied
Bernstein coefficients are metadata only. Conversion costs at most
`(d+1)(d+2)/2` rational multiply-adds per scalar coordinate (or use a precomputed
exact triangular transform for fixed `d`). The bound costs `d+1` absolute-value
comparisons and additions. Charge bit growth, transform arithmetic, and any
subdivision. This is not free cancellation: it only wins when the saved residual
integral or avoided search dominates this transform.

### Small strict-work example

Take one residual coordinate on `[0,1]`, `q(u)=1-u`, and target integral tolerance
`3/4`. The current coefficient-wise Taylor witness is
`|1|+|−1|/2 = 3/2`, so that candidate cannot pass. Bernstein coefficients are
`b_0=1,b_1=0`, giving `1/2`, so it passes. Suppose a frozen baseline pays 10
primitive operations to construct and check a replacement degree-one candidate
after failure, while Bernstein conversion and checking costs 6 operations. The
total is 6 versus 12 for the failed Taylor attempt plus replacement: a strict
saving with the same exact residual. This is a deliberately small accounting
counterexample; a real gate must compare against all supplied baselines, charge
failed candidates and transformations, and test coefficient scales where
Bernstein is worse. It does not establish a many-body advantage by itself.

## 2. A certified integrating factor on invariant weight sectors

Let `D(P)=-gamma*w(P)P` be the diagonal depolarizer and `C(P)=i[H,P]`. On a
subspace `V` spanned by Pauli words of one weight `w`, assume the checker verifies

```text
C(V) subseteq V.                                             (I1)
```

This is a finite graph check: for every Hamiltonian Pauli pair and every basis
word in `V`, recompute the commutator support and reject if a word leaves `V`.
On `V`, `D=-lambda I`, `lambda=gamma*w`, hence

```text
exp((C+D)t)|_V = exp(-lambda t) exp(Ct)|_V.
```

For Hermitian observables the Hamiltonian part is an isometry in operator norm,
so the residual contribution at the end of a segment is bounded by

```text
I_V = integral_0^h exp(-lambda*(h-u)) ||r_V(u)||_infty du.    (I2)
```

This is a valid replacement for the unweighted integral only after (I1) and
support of `r_V` in `V` have both been checked. Mixed-weight residuals are split
by sector and summed. If commutator support is not closed, use the unweighted
certificate; do not silently apply (I2). The same formula applies to every
intermediate endpoint, with `h` replaced by elapsed time.

For an exact scalar/anticommuting envelope
`||r_V(u)|| <= sum_k c_k u^k`, use

```text
I_V <= sum_k c_k J_k(lambda,h),
J_k = integral_0^h exp(-lambda*(h-u)) u^k du.                 (I3)
```

`J_0=(1-exp(-lambda h))/lambda` and
`J_k = h^k/lambda - (k/lambda)J_(k-1)` for `lambda>0`; at `lambda=0`,
`J_k=h^(k+1)/(k+1)`. Since exponentials are not rational, the checker stores
an outward rational interval for `exp(-lambda h)`, obtained by a fixed
alternating Taylor enclosure after range reduction (or a declared monotone
rational enclosure). It then propagates intervals through (I3), always rounding
up. The witness contains `lambda,h,c_k` and the interval for every `J_k`.

The extra work is sector-closure support checks plus `O(d)` scalar recurrences.
The saving can be large when `lambda h` is large: for a constant residual of size
`c`, the old term is `ch` whereas (I2) is `c(1-e^{-lambda h})/lambda <= c/lambda`.
At `gamma=0` there is no saving and the mechanism must report its charged
overhead. At `gamma=1/5` and `2`, this is precisely the kind of fixed-rate
headroom experiment the v7 gate requires; gamma cannot be tuned by the proposer.

## Cost and validity requirements

The two mechanisms compose only when each local witness is independently checked.
For each candidate report primitive observations, Hamiltonian-pair products,
support-closure checks, coefficient multiply-adds, rational bit maxima, grouping
comparisons, transform operations, exponential-enclosure operations, rejected
attempts, and final checker operations. Compare complete cost against conventional
sparse Taylor, projected BFS, residual-adaptive, Arnoldi, and anticommuting
grouping at each fixed gamma and tolerance. A lower residual number with higher
complete cost is not a win. A candidate that cannot prove sector closure or an
outward Bernstein/`J_k` enclosure must refuse.

These are algebraic tightening mechanisms, not novel identities: (B) is a basis
change plus convexity, while (I2) is variation of constants after a verified
scalar diagonal factorization. Neither supplies a learner, a causal
`m1 -> m2` transfer, or physical applicability beyond the declared unital-CP
certificate assumptions.
