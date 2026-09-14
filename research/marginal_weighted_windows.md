# Weighted windows: stronger energy bounds and an exact compatibility obstruction

Latest interval: [correlated hopping-filter upper](marginal_hopping_filter.md)
gives[-0.611636,-0.5581139722] at one million sites, with the upper rounded outward.

Follow-up: [quantitative projector extension and eight-site upper tiling](marginal_projector_extendibility.md)
now turn the compatibility obstruction into explicit Gram inequalities and
tighten the million-site interval to[-0.611636,-0.5294757424], with the displayed
upper rounded outward. The results below document the preceding window stage.

Exact nonuniform local windows now improve the certified energy interval for the uniform open Hubbard chain at U=4,t=1. At one million sites, the accepted six-site window gives

\[
-0.611636\le E_0/N\le -0.51542744,
\]

where the displayed upper endpoint is rounded outward. The exact rational interval has width approximately0.096208555314 per site, about47% below the preceding uniform-window width0.182174555314. This remains a coarse energy interval. Its local PSD matrices are at most200 by200, covering all4,096 local Fock states.

The four-site calculation also yields a sharp obstruction: its limiting physical state satisfies exact three-site overlap consistency, yet has no five-site extension that repeats the same four-site marginal. The accepted six-site positivity certificate supplies a linear energy separator against repeating its adjacent-pair data on a larger periodic chain.

## What the profile extension verifies

A local centered window now has the form

\[
K= -\sum_{i,\sigma}t_i(c^\dagger_{i\sigma}c_{i+1,\sigma}+\mathrm{h.c.})
+\sum_iU_i\left(D_i-\frac{n_i}{2}+\frac12\right).
\]

Profiles must be exact, bounded, nonnegative, reflection symmetric, and have their declared exact means. The CAR action, both reflection characters, Gram norms, and full local Fock coverage are still recomputed. The existing PSD matrix limit remains unchanged.

For a length-L overlap window representing the bulk parameters U,t, require

\[
\sum_iU_i=(L-1)U,\qquad \sum_it_i=(L-1)t.
\]

Summing all cyclic translates and dividing by L-1 then reproduces the centered periodic bulk operator exactly. Consequently a local all-Fock bound K≥aI gives the open-chain lower bound

\[
E_0(N)\ge \frac{Na}{L-1}-2t.
\]

The input mean for the onsite profile is therefore U(L-1)/L. The upper trial-state tiling can use a different window size. Tiled physical blocks must retain uniform bulk coefficients; the verifier rejects using nonuniform overlap profiles as such blocks.

A nonuniform centered shift does not vanish on every half-filled local vector. Its expectation is explicitly retained. A test with doubled edge sites and empty middle sites checks the distinction: its physical onsite energy is1 for the test profile, but its centered-window energy is6.

## Exact four-site optimization limit

The two-parameter family is

\[
U_i=[a,6-a,6-a,a],\quad
 t_i=[b,3-2b,b],\quad
0\le a\le6,\quad0\le b\le3/2.
\]

Numerical discovery required nine evaluations and proposed a≈0.5313730334,b≈0.75. Exact positivity accepts the rational point a=531373/1000000,b=3/4 with local lower endpoint

\[
l=-2.040424675.
\]

A separate integer-amplitude physical witness gives the exact constant expectation

\[
q=-\frac{913313562302334423671008}{447609546039589715350301}
\approx-2.040424674547836.
\]

Its expectation of both affine parameter derivatives is exactly zero after integer rounding. Thus, for every allowed a,b,

\[
\lambda_{\min}(K(a,b))\le\langle\psi|K(a,b)|\psi\rangle=q.
\]

Together with the accepted local lower point, this proves

\[
-2.040424675\le\max_{a,b}\lambda_{\min}(K(a,b))\le q.
\]

The exact bracket width is about4.521640796×10^-10. This is a certified limit of the lower-certificate family. The upper endpoint q limits the strength of this family; the full-chain variational upper energy is supplied separately by physical product states.

The verifier evaluates the affine expectation from three freshly constructed CAR operators. For general submitted witnesses, it maximizes the affine expression over the four rectangle corners. It never relies on optimizer stationarity. A test modifies an actual witness amplitude and verifies that the sharp requested upper limit is refused.

## A larger family has the same obstruction

The two three-site reduced density matrices of this four-site pure witness agree entry for entry: each has324 nonzero entries, and every difference is exactly zero. Therefore

\[
\langle A_{123}-A_{234}\rangle_\psi=0
\]

for every Hermitian even three-site operator A. Such differences telescope under cyclic translation.

The four-site profile variation itself has this form. Relative to the uniform local window, put d=a-3 and e=1-b, define Q_i=D_i-n_i/2+1/2 and T_i=Σσ(c†_{iσ}c_{i+1,σ}+h.c.), and take

\[
A=d(Q_1-Q_3)+e(T_1-T_2).
\]

Its left-minus-right embedding produces the onsite variation [d,-d,-d,d] and hopping-operator variation [e,-2e,e]. Hence the rational lower point belongs to the larger family consisting of the uniform window plus arbitrary three-site boundary differences. The same witness upper q applies throughout that larger family. Adding more such four-site overlap multipliers cannot overcome this ceiling.

Initial tests of zero-sum spin-exchange and density-density bond corrections also produced exactly zero expectations. Those tests found no violation; the complete reduced-matrix comparison explains why.

## Why matching overlaps still fails globally

The four-site witness density matrix is pure by construction. Its three-site reduced-state purity is

\[
\operatorname{tr}\rho_{123}^2
=\frac{127531781344489004873291493574281422558430209981}
{400708611411535170060004537416419717262281581202}
\approx0.3182656367<1.
\]

Suppose a five-site state had this pure witness as its first four-site marginal. Positivity forces the extension to factor as

\[
\rho_{12345}=|\psi\rangle\langle\psi|_{1234}\otimes\sigma_5.
\]

Its shifted four-site marginal is then ρ_234⊗σ_5. Its purity is at most trρ_234²<1, so it cannot equal the original pure four-site witness. This proves the absence of a five-site extension with identical neighboring four-site marginals, despite exact matching of the existing three-site overlaps.

This is an explicit missing extension constraint, rather than a numerical convergence problem. It establishes a concrete boundary between local consistency and global compatibility for this witness.

The six-site certificate also supplies a linear separator for its repeated adjacent-pair data. Their bulk energy per site is q/3≈-0.680141558183, and each site has exact mean filling1. The six-site positivity proof requires periodic-chain energy per site at least-0.611634 for N>6. The strict difference is approximately0.068507558183 per site. Therefore those pair data cannot be repeated as a stationary physical marginal on those periodic chains. They remain realizable on their original four-site block.

The final six-site replay computes this separation from fresh local positivity and family-witness checks. It does not trust a submitted moment or validation table.

## Six-site improvement

The accepted rational six-site window uses

\[
U_i=[3/20,83/25,653/100,653/100,83/25,3/20],
\]

\[
t_i=[5/12,5/4,5/3,5/4,5/12].
\]

Their sums are20 and5, as required for the same uniform bulk U=4,t=1. Exact all-sector positivity accepts a=-3.05817. Its largest reflected PSD matrix has dimension200. A separate bounded numerical refinement gives a slightly stronger candidate minimum about-3.0581689584; no family-optimality claim is made for six sites.

At N=1,000,000, retaining the same physical six-site product upper witness gives:

| Lower-certificate window | Lower energy per site | Exact interval width per site, approximately |
|---|---:|---:|
| Uniform six-site | -0.697602 | 0.182174555314 |
| Optimized four-site | -0.680143558333 | 0.164716113647 |
| Redistributed six-site | -0.611636 | 0.096208555314 |

The upper energy per site is approximately-0.515427444686 in each row. Finite open-chain boundary corrections are retained; N=12 and64 are also replayed, and the code takes the better of the tiling and overlap lower bounds.

These profiles are a stronger local representation of the same uniform bulk model. Molecular transfer, sharp accuracy at inexpensive block size, and general physical-marginal representability remain unresolved. The next constraint to pursue is extendibility beyond the locally consistent four-site block, alongside a stronger globally physical upper state.

## Validation and artifacts

- Full regression:420 tests passed in295.398seconds for the profile and affine-family implementation.
- After adding the full three-site consistency and five-site extension check, all16 relevant focused tests passed in3.756seconds. The full regression result is retained with that chronology.
- Separate standard-library replay recomputes the family bracket, the original nine global intervals, and three improved four-site-window intervals.
- A separate standard-library run recomputes the six-site window and three global intervals, then verifies the stationary pair-data separator.
- Current accepted input and replay artifacts are under `results/marginal_graded_hubbard8/weighted_window_family/`; numerical refinements are labeled separately.

```sh
OPENBLAS_NUM_THREADS=1 python results/marginal_graded_hubbard8/discovery/weighted_window_family.py
python -S results/marginal_graded_hubbard8/discovery/weighted_window_family.py --replay
python -S results/marginal_graded_hubbard8/discovery/six_weighted_window.py
python -S -m unittest tests.test_marginal_local_hubbard_block tests.test_marginal_local_energy_chain tests.test_marginal_window_family_bound -v
```

Production implementations are `marginal_local_hubbard_block.py`, `marginal_local_energy_chain.py`, and `marginal_window_family_bound.py` under `experiments/`. All lower proofs retain the existing local state/matrix budgets and recompute the original CAR action.
