# Exact global energy intervals from local Hubbard positivity

The open, half-filled Hubbard chain at U=4,t=1 now has independently replayed two-sided energy bounds whose local certificate cost is independent of total chain length at fixed block size. At one million sites, six-site local certificates give

\[
-0.697602\le E_0/N\le -0.515427444685994\ldots.
\]

The lower endpoint is exact as written; the upper endpoint is rounded from an explicit rational product-state witness. The width is approximately 0.1821745553 per site. This is a coarse global energy certificate, not a high-accuracy solution of the chain or a universal chemistry compiler.

The earlier moment-only result has therefore gained a separate global positivity construction. This construction uses established local-cluster reasoning; the advance here is the concrete implementation, exact all-sector checking, and measured certificate hierarchy in this workspace.

## Local operator and the sector issue

For a block of L sites define

\[
K_L(U,t)=H_L(U,t)-\frac U2(\widehat N_L-L).
\]

The original H_L has open nearest-neighbor hopping -t and onsite repulsion U. The centered onsite operator U(D_i-n_i/2+1/2) is zero on singly occupied states and U/2 on empty or doubly occupied states. On the global half-filled sector, K_N=H_N.

A lower certificate must establish K_L≥l_L I in **every local Fock sector**. Global half filling does not force each window or disjoint block to be half-filled. Nor does a grand-canonical minimum equal the half-filled minimum without an additional argument. The implementation avoids that assumption by checking every local particle and spin population.

`marginal_local_hubbard_block.py` constructs exact CAR actions on all 4^L local states for L=2,4,6. It groups states by (N_up,N_down), then splits each group into both characters of spatial reflection. The reflection includes the sign of the occupied-mode permutation. The actual CAR action is checked for Hermiticity, number-sector preservation, and exact commutation with reflection.

For each reflection orbit, the basis is either one fixed coordinate or the two vectors |s>±φ(s)|R(s)>. Their Gram norms are 1 or 2. Both characters, including fixed-point character signs, are included. The matrix checked is the congruence

\[
E^T K_L E-l_LE^TE.
\]

Exact denominator clearing and the existing fraction-free integer PSD checker verify it. The dimensions of the checked blocks sum to the full 4^L dimension:

| Sites | Local Fock dimension | Nonempty reflection blocks | Largest PSD matrix |
|---:|---:|---:|---:|
| 2 | 16 | 14 | 2 |
| 4 | 256 | 46 | 20 |
| 6 | 4,096 | 94 | 200 |

No PSD dimension limit was raised. The local work is still exponential in L.

The first draft incorrectly used determinant labels rather than amplitudes in an orbit Gram norm. Root review replaced that implementation before any accepted result was recorded. A focused test now accepts a lower endpoint below the exact two-site energy 2-2√2 and rejects one just above it; an artificially enlarged Gram norm would fail that check. Additional tests compare every local projected matrix entry with separately decoded original CAR actions and check full-basis orthogonality.

## Disjoint blocks give a physical upper and a global lower

Write N=qL+r, with even r and 0≤r<L. Every tiled block has an explicit integer-amplitude half-filled trial vector, whose exact Rayleigh quotient is u_L or u_r. A product of these vectors is globally half-filled. Hopping between blocks has zero expectation because every factor has definite charge. Thus

\[
E_0(N)\le qu_L+\mathbf1_{r>0}u_r.
\]

The norm of one hopping bond, summed over both spins, is at most 2t. There are q-1 removed bonds if r=0 and q if r>0. Local positivity gives

\[
E_0(N)\ge ql_L+\mathbf1_{r>0}l_r
-2t(q-1+\mathbf1_{r>0}).
\]

This is an operator inequality before restricting to global half filling. For N=qL, the interval width per site is at most

\[
\frac{u_L-l_L+2t}{L}-\frac{2t}{N}.
\]

This formula is conditional only on the actual local certificates and upper witnesses. It does not claim that arbitrary accuracy can be reached with small L.

## Overlapping windows strengthen the lower bound

Take an auxiliary periodic N-site chain with N>L and set

\[
v=U\frac{L-1}{L}.
\]

Each cyclic window uses the centered open-block operator K_L(v,t), including onsite interaction v. Every periodic bond appears in L-1 windows and every site appears in L windows. Exact counting gives

\[
\frac1{L-1}\sum_{j=1}^{N}K_L^{(j)}(v,t)=K_N^{\rm periodic}(U,t).
\]

Cyclic embeddings are embeddings of even CAR operators. Fermionic reordering signs are part of the operators, rather than being replaced by unsigned tensor permutations. Local positivity on all local Fock sectors survives each embedding. If K_L(v,t)≥a_LI, then

\[
E_0^{\rm open}(N)\ge \frac{N}{L-1}a_L-2t.
\]

The final -2t pays for removing the periodic wrapping bond. The code takes the larger of this lower bound and the disjoint-block lower bound. Independent sparse CAR tests check the full periodic operator identity for (N,L)=(4,2),(6,4), and a 16-dimensional exact PSD check covers both signs of the two-spin wrapping bond.

## Recorded exact certificates

Numerical eigensolvers propose endpoints and integer vectors; they do not certify any endpoint. The retained lower endpoints, all accepted by exact replay, are:

| L | Physical local U | Physical l_L | Overlap local v | Overlap a_L |
|---:|---:|---:|---:|---:|
| 2 | 4 | -0.829 | 2 | -1.237 |
| 4 | 4 | -1.954 | 3 | -2.348 |
| 6 | 4 | -3.093 | 10/3 | -3.488 |

The rational physical upper quotients are approximately -0.8284271247, -1.9531453087, and -3.0925653195. For a million-site chain, the resulting intervals are:

| Block size | Lower energy per site | Upper energy per site, rounded | Width per site, rounded |
|---:|---:|---:|---:|
| 2 | -1.237002 | -0.4142135624 | 0.8227884376 |
| 4 | -0.7826686667 | -0.4882863272 | 0.2943823395 |
| 6 | -0.697602 | -0.5154274447 | 0.1821745553 |

Nine global certificates cover N=12,64,1,000,000 for each block size. These include remainder blocks. The smaller-N bounds have explicit boundary corrections and need not have monotonically decreasing widths as N changes.

The input certificate is15,618bytes; the verbose replay record is348,308bytes and repeats local receipts across chain results. All six unique local certificates were recomputed in a fresh standard-library replay, which then composed the nine global intervals. The recorded batch took 43.202 seconds on CPU. Repeated local certificates are cached only within that fresh batch; no submitted receipt or persistent validation cache is trusted. A six-site local replay with a denominator-100,000 lower endpoint took 39.387 seconds; denominator-1,000 replay took 23.612 seconds. These are individual measurements, not a general performance law. The chosen coarser local endpoints add little to the much larger cluster-boundary uncertainty.

## What this resolves and what remains

There is now a concrete inner family of global positive-operator certificates, paired with globally physical trial states. It establishes coarse, two-sided accuracy without expanding the full chain. Its relationship to the marginal cone is explicit: the local PSD conditions imply globally valid linear energy inequalities, while the product states supply realizable marginals and upper witnesses.

Sharpness remains the main obstacle. Independent local positivity allows each overlapping window to minimize its energy without requiring all those local minimizing states to coexist in one global state. Stronger overlap consistency or a controlled coarse-graining representation is needed to close that gap cheaply. Increasing L alone eventually hits the exponential local calculation. The present verifier deliberately stops at L=6.

A bounded numerical follow-up redistributes local coefficients without changing their translated sums. For L=4, onsite values [a,6-a,6-a,a] sum to12 and hopping values [b,3-2b,b] sum to3. At the rational point a=1/2,b=3/4, an independent root calculation over all256 states gives a local numerical minimum -2.0405151147, or -0.6801717049 after division by3. The nearby a=53/100 has numerical density -0.6801416150. Repulsion is smaller at the edges and larger in the middle. These are numerical candidates only: the accepted verifier currently covers uniform local profiles. The diagnostic is retained in `local_energy_chain/weighted_window_diagnostic.json`; exact profile qualification and positivity are the next bounded implementation target.

Follow-up: [weighted-window certificates](marginal_weighted_windows.md) now implement and exactly verify those profiles. A sharper four-site point has a certified family-wide limit; a six-site profile reduces the million-site interval width to about0.09621 per site. The four-site limiting witness matches three-site overlaps exactly but cannot extend to five sites with repeating four-site marginals. The paragraph above records the earlier numerical discovery stage.

The construction also uses a uniform, finite-range, one-dimensional Hubbard model. All-sector local checks remove dependence on a ground-spin theorem, but do not establish molecular transfer, general N-representability, thermal observables, or synthesis predictions. The earlier sharp H8 interval is retained separately and is much narrower; the new larger-chain intervals are not equal-accuracy replacements for it.

For related rigorous lower-bound hierarchies that retain consistency between coarse-grained reduced states, see [Kull, Schuch, Dive, and Navascués, Phys. Rev. X 14, 021008 (2024)](https://arxiv.org/abs/2212.03014). Their construction supports the choice of next target; it is not a validation of the implementation here.

## Reproduction and artifacts

```sh
OPENBLAS_NUM_THREADS=1 python results/marginal_graded_hubbard8/discovery/local_energy_chain.py
python -S results/marginal_graded_hubbard8/discovery/local_energy_chain.py --replay
python -S -m unittest tests.test_marginal_local_hubbard_block tests.test_marginal_local_energy_chain -v
```

- Local engine: `experiments/marginal_local_hubbard_block.py`.
- Global composition and fresh-batch replay: `experiments/marginal_local_energy_chain.py`.
- Proposals: `results/marginal_graded_hubbard8/local_energy_chain/certificate.json`.
- Accepted independent replay: `results/marginal_graded_hubbard8/local_energy_chain/independent_replay.json`.
- Focused validation: ten tests passed in 1.708 seconds; `local_energy_chain_focused_validation.log`.
- Full regression: **414 tests passed in292.684seconds**; `local_energy_chain_full_validation.log`.
