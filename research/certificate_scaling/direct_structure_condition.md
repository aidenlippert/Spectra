# Direct discovery by signed diagonal dominance

## Goal

The local-overlap theorem bounds a certificate once one is supplied. This
condition addresses the missing direction: a verifier can inspect a sparse
operator dictionary and construct an SOS certificate when its coefficient
matrix is signed diagonally dominant. It is a screening rule, not a claim
about generic molecular Hamiltonians.

## Input condition and recognizer boundary

Fix number-conserving fermionic words `v_1,...,v_s`, each with a recorded CAR
normal form and support of at most `w` modes. Let `P_N` be the fixed-particle
projector. The input supplies a rational identity

```
P_N (H - b I - J - R) P_N = P_N (v^dagger G v) P_N,
```

where `v^dagger G v = sum_(i,j) G_ij v_i^dagger v_j`, `G` is real symmetric,
and `P_N J P_N=0` is checked explicitly. The residual is `R=sum_X R_X` with
rational local norm bounds `u_X >= ||P_N R_X P_N||`.

The directly checkable condition is

```
G_ii >= sum_(j != i) |G_ij|       for every i,       G_ii >= 0.
```

The recognizer is efficient only when `H` is already supplied in this sparse
word dictionary (or a CAR parser can map its terms there). It builds `G` by
coefficient matching, then checks the inequalities. Thus `G` is not a magical
certificate input: the structural hypothesis is that the Hamiltonian's chosen
operator coordinates expose a sparse diagonally dominant Gram block. Generic
molecular two-electron tensors need not have this form.

## Constructive theorem

If the condition holds, construct the exact factor-width-two SOS

```
v^dagger G v = sum_(i<j) |G_ij| (v_i + sign(G_ij) v_j)^dagger
                              (v_i + sign(G_ij) v_j)
                 + sum_i d_i v_i^dagger v_i,
d_i = G_ii - sum_(j != i) |G_ij| >= 0.
```

Consequently, `E_0(H;N) >= b - sum_X u_X`. Positivity survives fixed-sector
restriction and the only discarded equality is the explicitly checked `J`.

### Proof

Expanding one pair term contributes `|G_ij|` to both diagonal entries and
`G_ij` to each symmetric off-diagonal entry. The diagonal leftovers are
exactly `d_i`; summing gives `G` entry by entry. Every summand is positive, so
its fixed-`N` expectation is nonnegative. The residual triangle inequality
contributes at most `sum_X u_X`.

## Exact construction cost

For `e` nonzero off-diagonal entries, construction emits at most `e+s` SOS
atoms. Dominance costs `O(s^2)` rational comparisons for a dense matrix and
`O(e+s)` for sparse adjacency lists. Expanding each two-word atom and checking
the CAR map costs `O(e*poly(w))` word operations, or `O(e*t)` with a sparse
normal-form table having at most `t` output monomials per product.

The residual cost is separate. If each residual block has at most `r` local
modes, an exact local norm bound may cost `2^(O(r))` by local matrix
construction or interval eigenvalue bounding. The theorem does not hide this
factor. With `q` residual supports, the additive error is `eta=sum_X u_X` and
the total verification cost is

```
O(e*sparse_map_cost + q*2^(O(r)) + bitcost of all rational data).
```

## Recognizable interacting family

Take pair-annihilation words `v_(ij)=a_j a_i` on the edges of a sparse graph
and let `G` be a rational signed diagonally dominant matrix on those edges.
Then `sum G_(ij),(kl) v_(ij)^dagger v_(kl)` is a quartic, number-conserving
Hamiltonian containing pair-density and pair-exchange terms. A bounded-degree
graph gives `e=O(s)` for sparse `G`; the certificate is constructed from the
Hamiltonian coefficients in one pass over its edges. This family is physically
recognizable but excludes arbitrary long-range Coulomb tensors and many
strongly correlated active spaces.

## Decisive failure boundary

The rule is sufficient, not necessary. For

```
G = [[1, 2], [2, 4]],
```

`G` is positive semidefinite and has the one-atom SOS
`(v_1+2v_2)^dagger(v_1+2v_2)`, but row one fails dominance. Thus the rule can
reject compact certificates found by another Gram factorization. If `G` is
indefinite, no PSD Gram certificate in this fixed dictionary can represent it
without a residual or additional words. The condition is basis dependent and
is a constructive pricing filter, not a universal fermionic RIP.

## Minimal exact experiment

Use `v_1=a_2a_1`, `v_2=a_4a_3` and `G=[[3,-1],[-1,2]]`. Dominance holds, with
`d=(2,1)`, so construction emits

```
(v_1-v_2)^dagger(v_1-v_2) + 2 v_1^dagger v_1 + v_2^dagger v_2.
```

An exact stdlib checker can compare the rational coefficient entries after
expansion and report `eta` from separately supplied residual blocks. This
checks construction in coefficient space; a global Fock matrix is unnecessary.

The construction uses rational weights multiplying unnormalized two-word
squares. It does not require taking square roots of `|G_ij|`, which would
generally leave the rational coefficient field.
