# Quantitative extendibility breaks the six-site window ceiling

An exact overlapping-projector constraint and a retuned local operator
improve the lower energy bound for the half-filled open Hubbard chain
with U=4, t=1, V=1/2. At one million sites the final bound is

\[
E_0/N\ge-63102110090479/97656250000000
\simeq-0.646165607326505.
\]

The preceding lower bound was -0.6632065. More significantly, the new
value exceeds -0.663023831373854, the exact ceiling previously proved
for every six-site profile with the required coefficient sums, including
arbitrary even five-site telescoping boundary corrections. Thus this
improvement cannot be obtained by further tuning that old family.

A fresh physical transfer-state replay gives the upper bound. The final
ground-energy interval, with displayed endpoints rounded outward, is

\[
\boxed{-0.64616560733\le E_0/N\le-0.61067634705.}
\]

Its width is approximately 0.0354892602753 per site, closing 32.4402% of
the previous interval. These are energy units with t=1. The result applies
to the specified model and chain; it is not an error guarantee for chemistry.

## The physical information added

Let phi be the physical six-site cyclic state with 400 integer amplitudes
from `density_transfer/window_ceiling_certificate.json`, and let
P=|phi><phi|/<phi|phi>. Its five-site reductions agree exactly but are
mixed. A global state cannot have this same pure six-site marginal on
both adjacent windows. The first accepted certificate quantified this obstruction:

\[
P_0+P_1+P_2+P_3\preceq B I,
\qquad B=1084870113/500000000=2.169740226.
\]

This inequality holds on the full nine-site Fock space. It is verified
without constructing a 262,144-dimensional many-body operator. Stack the
isometries that insert phi into each window, normalized by its norm. Their
Gram matrix has the same nonzero spectrum as the projector sum. Its 256
columns split by actual conserved spin numbers into matrices of dimension
at most 36. Exact integer positive-semidefinite checks certify the ceiling.
Every environment configuration is included.

On a periodic N-site chain, summing all translated copies counts each P_i
four times. Consequently

\[
\frac1N\sum_iP_i\preceq\theta I,
\qquad\theta=B/4=0.5424350565.
\]

Cyclic copies are conjugated by the physical fermionic mode-translation
unitary, retaining signs across the closing bond. The projector is parity
even, so contiguous embeddings use the ordinary local tensor matrix.

## Conversion into a ground-energy certificate

The centered six-site local operator K has reflected profiles

```
onsite:  [3/20, 83/25, 653/100, 653/100, 83/25, 3/20]
hopping: [5/12, 5/4, 5/3, 5/4, 5/12]
density: [5/24, 5/8, 5/6, 5/8, 5/24]
```

Their sums are 20, 5, and 5/2. The accepting local inequality is

\[
K+\kappa P\succeq\ell I,
\qquad\kappa=30527/200000,
\quad\ell=-31634217/10000000.
\]

It is checked in all 94 reflected spin-number sectors, covering all 4,096
local Fock states; the largest matrix has dimension 200. Reflection-orbit
columns are not normalized. The verifier uses their exact Gram factors
when subtracting ell*I and forms the projector by exact congruence.

Translated K windows sum to five times the periodic target Hamiltonian,
up to the charge-centering term that vanishes at global half filling.
Combining the two operator inequalities yields

\[
e_{\rm periodic}\ge\frac{\ell-\kappa\theta}{5}
=-0.6492432549697755.
\]

Removing the closing hopping and density bond costs at most
(2t+|V|)/N=0.0000025, giving the first open-chain bound
-0.6492457549697755. This alone closed 26.5766% of the old interval.

## The stronger five-projector certificate

Five adjacent projectors give the exactly accepted inequality

\[
\sum_{j=0}^{4}P_j\preceq\frac{2621073531}{10^9}I,
\qquad\theta=0.5242147062.
\]

The support is ten sites. The full Gram has 1,280 columns and largest
spin block dimension 180. Its dimension counts independently satisfy
5*binomial(4,Nup-3)*binomial(4,Ndown-3) in each spin sector, including
both exterior environments for every window position.

A bounded 118-evaluation profile search proposed a new local operator.
Fresh exact CAR matrices agreed with its affine numerical model before
the final positivity replay. With profile notation

```
onsite:  [a, b, 10-a-b, 10-a-b, b, a]
hopping: [p, q, 5-2p-2q, q, p]
density: [d, e, 5/2-2d-2e, e, d]
```

the final parameters are

| a | b | p | q | d | e |
|---:|---:|---:|---:|---:|---:|
| 0.204105 | 3.039489 | 0.478826 | 1.260412 | 0.147254 | 0.372037 |

The exact penalty is kappa=3336/15625=0.213504 and the accepted local
lower value is ell=-3898617/1250000=-3.1188936. All 94 local sectors
again pass. Substitution into (ell-kappa*theta)/5-(5/2)/N gives the final
lower endpoint above. Search optimality is not asserted.

The upper replay reconstructs the physical eight-site source and evaluates
the same target using the gate F=I-0.374704*h+0.177582*h^2. An exact 24-site
contraction is contained in the independently rounded transfer enclosure.
The million-site transfer uses 160-bit outward arithmetic; its trial-energy
enclosure width is below 1.3e-37 per site. The much larger ground-energy
interval is uncertainty in the variational and positivity bounds.

## Exact arithmetic and verification

Two exact local replays using different denominator treatments produced
identical complete receipts. The direct denominator-scaled version took
299.349 seconds; the rank-one-aware version took 96.485 seconds in these
runs. These are observed runtimes, not a controlled scaling benchmark.

For a rational base matrix A and rank-one addition r*ww^T/s, the optimized
check retains the positive initial divisor s in fraction-free elimination.
Every update must divide exactly. Each positive-pivot update is a positive
multiple of a Schur complement; zero pivots require a zero coupling row.
This avoids carrying unnecessary denominator powers while preserving the
PSD criterion. Tests cover singular, zero-row and indefinite cases, invalid
divisors, and nonexact-division refusals.

The first focused suite passed 14 tests; the final focused suite passes
15 tests in 4.956 seconds. It includes an independent two-window
tensor embedding of the actual cyclic source, all Gram dimensions, and
target, spin, reflection, penalty and ceiling refusals. The initial full
regression passed 468 tests in 336.785 seconds. The final regression passed
all 469 tests in 347.596 seconds. Logs and receipts are listed in
`results/marginal_final_validation.json`.

The four-projector matched replay completed in 101.592 seconds. The refined
five-projector replay completed in 290.226 seconds. All 38 replay source/input
hashes were checked. A subsequent correction of the driver's descriptive
support-size label is documented with its original source snapshot and an
AST comparison confirming that only the descriptive string changed.

Two numerical agent outputs were rejected: one omitted the right exterior
environment when assigning spin sectors; another shifted affine derivative
indices and exceeded its evaluation cap. Both are explicitly preserved as
rejected diagnostics. The corrected profile probe enforces a hard cap,
retains reflection Gram normalization, and reconstructs the final physical
matrix independently of the affine interpolation. Only exact accepting
replays contribute to the reported energy bounds.

```sh
OPENBLAS_NUM_THREADS=1 python -S results/marginal_graded_hubbard8/discovery/six_site_projector.py --certificate refined_certificate.json
```

## What remains open

This supplies an explicit separating inequality missing from local overlap
consistency. It does not characterize the full cone of physical marginals.
We still need a method that discovers enough such constraints at controlled
cost, closes the certified energy interval to a requested tolerance, and
transfers beyond these specified one-dimensional interactions. Thermal
states, dynamics, response and synthesis require further constructions.

**Continuation update:** the separate-penalty proposal below is historical and
structurally dominated. See [the joint-projector analysis](marginal_joint_projectors.md)
for the proof and the current coupled attack.

The refined numerical local minimum occurs in the charged (Nup,Ndown,
reflection)=(2,3,+1) sector. A concrete next family adds projectors for this
five-particle state and its spin-swapped and particle-hole-related partners,
while retaining the existing half-filled penalty. Separate exact overlap
ceilings would allow the lower formula

\[
e_{\rm periodic}\ge
(\ell-\kappa\theta_{\rm half}-\lambda\theta_{\rm charged})/5.
\]

For four charged source sectors (2,3),(3,2),(4,3),(3,4) and four windows,
the Gram would have 1,024 columns and maximum spin-block dimension 96.
Each local spin sector would still contain at most one rank-one penalty,
preserving the useful denominator structure. These counts are exact
combinatorics; no charged-state energy improvement has yet been computed
or certified. Numerical proposals remain separate from exact certificates.

For a fixed odd-particle source, a graded embedding changes each whole
environment column by the constant sign (-1)^(Nleft*Nsource). This gives
Vgraded=Vunsigned*D for diagonal signs D, hence identical VV-dagger and
unitarily congruent Gram matrices. Odd source parity therefore does not
invalidate the unsigned projector-sum embedding. An initial contrary
audit claim was corrected. Actual source symmetry, local positivity and
the new charged overlap ceiling still require verification.
