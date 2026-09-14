# Local cluster certificates: what a logarithmic radius would require

## The tempting claim

For a one-dimensional finite-range chain with a unique gap `Delta`, it is
tempting to assert that a block certificate of radius `R=O(log(M/epsilon))`
gives an absolute ground-energy error at most `epsilon`. The gap and exponential
decay of local correlations do not prove this. Ground energy is extensive: an
error `q_R` per translated block gives total error of order `M q_R`.

## A rigorous conditional route

Consider the perturbative family

```
H(lambda) = H_parent + lambda V,
H_parent = sum_x h_x,
```

where `h_x` is a two-mode dimer Hamiltonian with a unique product ground state,
local gap at least `Delta`, and `V` is nearest-neighbor, number-conserving with
`||V_e|| <= v`. Suppose a linked-cluster construction supplies an explicitly
verified bound

```
|e_0(lambda)-e_0^[R](lambda)| <= C*rho^(R+1), rho<1,
```

for the thermodynamic energy density, plus a finite-size boundary bound
`C_b*rho^R`. Here `e_0^[R]` is assembled from connected clusters of radius at
most `R`; every cluster is replayed by a local exact SOS calculation. Then

```
|E_0^(M)-[M*e_0^[R]+boundary_R]|
    <= M*C*rho^(R+1)+C_b*rho^R.
```

Choosing `R >= log((M*C+C_b)/epsilon)/|log(rho)| - 1` gives absolute error at
most `epsilon`. This is a sufficient theorem, but its hard hypothesis is the
certified linked-cluster remainder. A spectral gap alone does not provide it.

For a perturbative proof one could derive `rho` from a convergent polymer
criterion involving `q*v/Delta`, with `q` the interaction degree. The criterion
must include complex-energy denominators and polymer combinatorics. Writing
`rho=lambda*v/Delta` without proving those counts is not a certificate.

## Cost and compatibility

At radius `R`, a connected cluster has `O(R)` dimers and `O(R)` boundary
interfaces. A local exact calculation costs roughly `2^(O(wR))`; the number of
translated cluster types is polynomial in `R` for a translation-invariant
chain. Boundary redistribution requires an exact overlap map: contributions
must agree on shared subclusters and signed inclusion-exclusion coefficients
must be replayed exactly. A local lower bound without this compatibility map
cannot be summed into a global lower bound.

The total bound is `M*C*rho^(R+1)+C_b*rho^R`, not one independent `epsilon` per
cluster. With direct residuals the safe bound remains `sum_X u_X`, which is
extensive at fixed interaction strength.

## Why the gap alone is insufficient

Even a gapped chain can have an extensive first-order energy response. For a
product parent ground state and `V=sum_e V_e`,
`dE_0/dlambda|_0=<Psi_parent|V|Psi_parent>` is generally `Theta(M)`. A gap
controls stability and often local correlations; it does not make omitted
energy contributions small. The cluster certificate must compute extensive
low-order terms exactly before bounding a decaying remainder.

Moving half of a boundary interaction to each adjacent block does not make block
inequalities compatible unless the shared boundary is included in an exact
overlap SDP or paid for by its norm.

## Concrete frontier experiment

Use the dimer parent with nearest-neighbor pair hopping from
`cluster_structural_check.py`. Construct connected cluster energies through
radius `R`, normalize by exact inclusion-exclusion, and record a rational
interval for `rho_R=|I_(R+1)|/|I_R|`. The route succeeds only if an interval
proof establishes `sup_R rho_R<1` and a separate boundary enclosure. Otherwise
it is an empirical decay plot, not a global certificate. The fixed-strength
norm-tail check already shows why raw edge norms cannot solve this problem.

## Scope

The parent/product case, commuting chains, and free-fermion diagonalization are
controls, not evidence for a new interacting theorem. The difficult step is a
checkable linked-cluster remainder with compatible boundary redistribution.
Without it, bounded local certificates plus a finite gap do not imply an
absolute chemical-accuracy interval at logarithmic radius.
