# Exact-solvable reference expansion: bounded preflight

## Algebra and certificate

Let `G=G0+B` on the finite-dimensional Hermitian operator space. Assume both
`U(t)=exp(tG)` and `U0(t)=exp(tG0)` are positive, unital, and operator-norm
contractive, and set `D0(t)=U0(t)O`. Define recursively

```text
Dk(t) = ∫_0^t U0(t-s) B D{k-1}(s) ds,
Pm(t) = Σ_{k=0}^m Dk(t).
```

Differentiation gives `D0'=G0D0` and
`Dk'=G0Dk+B D{k-1}`. Hence

```text
Pm' - G Pm = -B Dm,
```

exactly. The sign in the proposed lemma is therefore correct. If
`K >= ||B||_{∞→∞}`, contraction implies inductively
`||Dk(t)||∞ <= ||O||∞ (K t)^k/k!`. Applying the same residual/variation of
constants certificate as v7 gives

```text
||U(T)O-Pm(T)||∞
 <= ∫_0^T ||B Dm(s)||∞ ds
 <= ||O||∞ (K T)^(m+1)/(m+1)!.
```

There is no extra `exp(KT)` factor: that factor would arise from bounding the
perturbed propagator by a generic growth estimate, whereas the hypothesis here
already gives contraction of the full `U(t)`. A looser estimate with
`exp(KT)` remains valid but is unnecessary.

This is a finite-dimensional bounded-perturbation/variation-of-constants
construction; no priority claim is made. A verified primary overview of
exponential-integrator construction and analysis is
[Hochbruck and Ostermann, 2010, author-hosted PDF](https://na.math.kit.edu/download/papers/acta-final.pdf).
The specific factorial estimate above is proved directly here from the two
contraction assumptions. An agent-supplied Desch–Schappacher citation could
not be independently verified in this pass and is not used as support.

## Exactly solvable reference and representation

For `G0` consisting of one selected Pauli-axis field per qubit plus uniform
local depolarizing terms, local Clifford rotations map each selected axis to Z.
The one-qubit eigenoperators are `{I,Z,σ+,σ−}`, with exponents
`0`, `−gamma`, `−gamma+2ih`, and `−gamma−2ih`. Product exponents add.
An arbitrary direction hx X+hy Y+hz Z cannot generally be mapped to Z by a
Clifford operation; it would require a general rotation and extra algebraic
number accounting. Other axis components can instead remain in B. Thus `U0(t)` acts diagonally by exact exponential
factors, and repeated Duhamel integration produces exponential-polynomial
terms. Coincident frequencies must be combined with the correct limiting
polynomial factors; treating equal frequencies as distinct denominators is an
algebraic failure mode.

An independent checker could accept a finite sum of such terms with exact
rational/complex-rational frequencies and coefficients, while evaluating the
real exponential/trigonometric factors with outward bounds. This is a design
option only; the trusted v7 checker remains unchanged in this preflight.

For `B=i[V,·]` with Hermitian `V`,
`||B(X)||∞ = ||i[V,X]||∞ <= 2||V||∞||X||∞`, so `K=2||V||∞` is valid. The
bound is generally conservative and must not be replaced by a scalar spectrum
of a nonnormal generator.

## Resource accounting and decision boundary

At order `k`, each Duhamel application requires applying `B` to the preceding
term and integrating against the known `G0` modes. Cost is governed by the
number of retained operator modes and by collisions that create polynomial
degree, plus exact coefficient bit lengths and outward exponential bounds.
The candidate and full dynamics share the same `H`, `γ`, and horizon `T`; all
mode generation, convolution, repeated-frequency handling, failed orders, and
certificate replay must be charged. The residual certificate can be evaluated
with the existing Pauli majorants after conversion, but a new exponential-term
checker would itself require a separately audited implementation.

The factorial tail is a rigorous baseline, not a headroom claim. For the
current V7 family, G0 removes onsite Z fields and depolarization, but the
remaining interactions are not automatically small. For n6 XXZ, the simple
triangle bound is K<=8.8; with T=.2 this gives KT<=1.76. In mixed cases the
other onsite axis and YY interactions also remain. No cost win follows from
calling the reference exactly solvable. No claim of
solved learning, novelty, or superiority follows from this lemma. A preflight
should compare accepted certified error and total exact work against v7 Taylor
and Arnoldi under identical `H,γ,T`, and refuse when mode growth, coefficient
growth, or outward exponential bounds erase the apparent tail advantage.


## Root review of the proof

Each complete D_k is Hermitian when O is Hermitian and both generators preserve adjoints, although individual exponential modes may be complex. The contraction proof applies to these complete Hermitian sums. A checker must verify conjugate pairing or otherwise establish Hermiticity; it cannot apply that premise independently to an arbitrary unpaired complex coefficient. All physical assumptions and the final error target must remain identical across comparison arms.
