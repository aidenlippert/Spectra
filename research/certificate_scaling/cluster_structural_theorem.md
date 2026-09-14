# A compositional cluster theorem with an interacting residual budget

## Setup

Partition `M` fermionic modes into disjoint clusters `C_1,...,C_p`, each of
size at most `w`. Let `n_x` be the particle number in cluster `x` and write

```
H = sum_x H_x + sum_(e in E) V_e,
```

where `H_x` acts only on `C_x`, and `V_e` acts on the union of clusters in
edge `e`, with `|e| <= r`. Assume rational local certificates, for every local
particle number `n` that can occur,

```
P_(x,n) (H_x - b_x(n) I) P_(x,n) >= -eta_x(n) I.
```

The local bounds may be verified by a compact SOS identity and a local residual
norm computation. No global Fock basis is needed. For each interaction edge,
provide a rational bound

```
||P_N V_e P_N|| <= u_e.
```

The unrestricted local norm is safe but may be loose; a fixed-total-`N` bound
is allowed only when its projection is explicitly justified.

## Theorem (cluster composition)

For fixed total particle number `N`, define

```
B_N = min_{n_1+...+n_p=N} sum_x b_x(n_x),
Eta_N = max_{n_1+...+n_p=N} sum_x eta_x(n_x),
U = sum_e u_e.
```

Then

```
E_0(H;N) >= B_N - Eta_N - U.
```

If all local residuals are uniform, `eta_x(n)<=eta_x`, then one may replace
`Eta_N` by `sum_x eta_x`. The interaction budget is the *sum* of edge norms;
assigning an `epsilon` allowance independently per edge without summing it is
not a valid total-error statement.

### Proof

The fixed-total sector decomposes as the orthogonal direct sum over local
occupancies `(n_x)` with sum `N`. On each block, the local certificates add,
giving a lower bound `sum_x b_x(n_x)-sum_x eta_x(n_x)`. The minimum over blocks
is at least `B_N-Eta_N`. The interaction sum obeys

```
<psi|sum_e V_e|psi> >= -sum_e ||P_N V_e P_N|| >= -U.
```

Taking the infimum proves the claim. This is a triangle-inequality theorem;
overlapping interactions need not commute.

## Exact disjoint baseline and upper witness

When `E` is empty, the lower bound is the exact product-cluster baseline if
local certificates are exact (`eta_x(n)=0`):

```
E_0(sum_x H_x;N) = min_{sum n_x=N} sum_x E_0(H_x;n_x).
```

An upper bound for the interacting system can be made without global
enumeration by choosing local normalized witnesses `phi_x` with occupancies
`n_x` and forming their graded tensor product `Phi`. Its energy is evaluated as

```
<Phi|H|Phi> = sum_x <phi_x|H_x|phi_x> + sum_e <Phi|V_e|Phi>,
```

using only edge-local reduced data. Thus an energy interval can be reported
with total width `U + Eta_N + (upper-witness error)`, provided every local and
edge expectation is exact or interval certified.

## Scaling and precision

Let `q` be the maximum edge degree, `s_x` the local certificate dictionary
size, and `L` the maximum rational bit length. Local verification costs at most
`p * poly(s_x) * 2^(O(w))` using exact local matrices; each edge residual costs
`2^(O(wr))` in the worst case. Computing `B_N` by dynamic programming over the
integer occupancies costs `O(p N^2)` for bounded cluster occupancy (or the
corresponding convolution cost), without constructing `binom(M,N)` states.
The printed additive error is a rational sum; its bit length grows linearly in
the input bit lengths plus `O(log p)` additions. If `|E|=O(pq)`, the coupling
budget is extensive unless `u_e` decreases with system size.

The theorem is therefore useful for fixed local complexity and weak aggregate
coupling. It does not claim a size-independent global chemical-accuracy bound
for an extensive interacting material.

## Nontrivial finite family

Take clusters to be two-mode dimers with local pair penalties and let adjacent
dimers have pair-exchange interactions `V_e` of rational amplitude `g_e`.
Each dimer has an exact local SOS certificate. For a two-dimer edge, the pair
exchange operator has norm at most `|g_e|` (or a directly verified small-block
norm), so the chain lower bound is the exact dimer occupancy baseline minus
`sum_e |g_e|`. A product of local dimer witnesses supplies the upper endpoint,
with each exchange expectation evaluated on two dimers only. This is genuinely
interacting for nonzero `g_e`, while its certificate construction remains local.

## Falsifiable discovery test

For localized active-space calculations, measure the tuple

```
(max cluster certificate width, max edge support, sum_e u_e, observed interval width)
```

as cluster size and chain length grow. The structural hypothesis predicts that
local widths remain bounded and that the aggregate edge budget is the dominant
error. It fails if local residuals or edge norms grow with cluster size, if
occupancy-sector certificates cannot be reused, or if the sum of edge budgets
remains too large for the target accuracy. A failure is informative: it rules
out this compositional route for that family without saying that compact
nonlocal certificates are impossible.
