# Exact symmetry-orbit moments and a sharper H8 energy certificate

The original eight-site, half-filled open Hubbard chain at U=4 and t=1 now has two independently replayed ground-energy intervals:

| Complement discovery | Exact lower endpoint | Conservative upper endpoint | Exact interval width, approximately |
|---|---:|---:|---:|
| Fresh range-two charge-product discovery, gamma=-3.81 | -4.23585 | -4.23580593 | 0.00004406075472 |
| Enumerated range-three discovery, gamma=-3.76 | -4.235845 | -4.23580593 | 0.00003906075472 |

The common exact upper endpoint is

\[
-\frac{120558410244999997006043477431315508062486123707}
{28461740687412294825487447957163629238443780483}.
\]

These are hopping-energy units for the specified lattice model. They are not molecular benchmarks or material-property predictions. The ground-spin step retains the explicit dependency on [Lieb's repulsive half-filled bipartite Hubbard theorem](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.1201), with the original model hypotheses checked.

## What changed

The preceding energy replay expanded 14 response vectors into thousands of determinant amplitudes. The new replay recomputes full moments M_n=VᵀHⁿV in signed symmetry orbits. It then uses the exact nonorthogonal renewal recurrence to recover the projected QHQ moments. The lower bound needs only small rational matrix operations after that contraction.

The response is one degree-ten polynomial applied to all 14 leakage vectors. The upper witness is a degree-eleven polynomial combination of the 14 valence-bond vectors, represented by 168 integer coefficients in the fixed Chebyshev basis T_k((H-8I)/12). That change avoids rounding a long, expanded determinant vector. Its norm and energy are recomputed exactly from the same full moments.

The strongest calculation uses M_0 through M_24. These require powers only through H^12 because Hermiticity gives

\[
M_n=(H^{\lfloor n/2\rfloor}V)^\top H^{\lceil n/2\rceil}V.
\]

The standalone moment API was already bounded at order 24. The new moment-energy family admits response degree at most ten and upper degree at most eleven within that bound. The older determinant-vector energy verifier keeps its original polynomial limit. No determinant-action cap was raised.

In the recorded run, all 25 moments took about 1.44 seconds. Full sharp energy replay took about 3.52 seconds with the range-two gap and 3.74 seconds with range three. These are individual CPU-run measurements, not equal-work speedup estimates against previous code.

## Exact signed-orbit construction

Spin exchange F swaps the two spin modes at each site and contributes the CAR phase (-1)^D. Particle–hole conjugation C sends each creation operator to (-1)^site times its annihilation operator. On a determinant s its phase is

\[
(-1)^{\sum_{i\in\mathrm{occupied}(s)}(\lfloor i/2\rfloor+i)}.
\]

The implementation checks spin-exchange invariance symbolically. It also checks the exact particle–hole identity modulo fixed particle number. For this model,

\[
CHC^{-1}-H=4(8-\hat N),
\]

which vanishes on the specified half-filled sector. This is not advertised as equality on the entire Fock space. The module additionally requires balanced spin populations and checks the invariance of each supplied boundary vector. A negative stabilizer makes an orbit forbidden.

For an allowed orbit, use the unnormalized vector E_r=Σ_t φ_r(t)|t⟩. Its squared norm is the orbit size. Applying H to one representative gives the quotient coefficient

\[
\widetilde H_{ur}
=\frac{|O_r|}{|O_u|}\sum_{t\in O_u}\phi_u(t)H_{tr}.
\]

Inner products carry the corresponding orbit-size weight. Omitting this normalization factor would invalidate the contraction.

The energy verifier never accepts a submitted moment table. It recomputes the Hamiltonian actions, moments, renewed Q moments, positive response solve, Schur pivots, and physical polynomial Rayleigh quotient. It also reruns the complete 14-dimensional singlet embedding and the charge-DP complement certificate.

## Verification and failed alternatives

- Every signed-orbit action on the four-site chain matches independently expanded CAR actions. Orbit-size-weighted Hermiticity is checked.
- Four-site full moments through order 14 match the direct route exactly.
- Every one of the 2,940 H8 matrix entries in M_0 through M_14 matches ordinary determinant CAR contraction, using 3,870 original source actions versus 995 quotient source actions at that order.
- A nonuniform four-site bipartite ring with an additional hopping of 2/3 passes exact moment comparisons through order eight.
- The nonorthogonal renewal identity is independently checked through projected moment order 12 on four sites.
- Tests reject broken spin or particle–hole symmetry, non-invariant boundaries, floating boundary coefficients, invalid upper witnesses, unsupported moment orders, dependent response blocks, and unproved lower endpoints.

A targeted local optimization of the degree-seven response failed to certify -4.236. Its exact rejection is retained. Adding a second response block, for 28 directions total, also did not reach that endpoint in the tested proposals. These are finite unsuccessful searches, not impossibility proofs for those families. Higher degree within the new moment representation supplied the accepted sharper interval.

The range-three count-profile proposal has 91 parameters and numerical optimum approximately -3.74936844. Exact DP accepts gamma=-3.76. Discovery enumerated 1,106 charge patterns, and the DP's largest local window contains all eight sites. The fresh range-two alternative is retained separately to keep those distinctions visible.

## The remaining scaling obstruction

At moment order 24, H8 uses 1,239 distinct representative determinant actions, maximum vector support of 1,250 orbits, and an orbit-cache footprint covering 4,864 determinant labels. The ordinary CAR images from representative actions reference 2,169 labels; the larger orbit-cache number is also reported rather than hidden. The full simultaneous positive-character symmetry sector has dimension (4900+70+70)/4=1260 by the finite-group trace formula, so the high-order vectors still occupy almost that entire symmetry sector.

The quotient therefore removes a finite symmetry factor. It does not establish polynomial scaling.

Single product-of-singlets boundary probes make the limit concrete:

| Sites | Moment order | Result | Representative source actions |
|---:|---:|---|---:|
| 8 | 24 | Accepted | 1,238 |
| 12 | 4 | Accepted | 208 |
| 12 | 8 | Refused at existing action cap | 4,096 |
| 16 | 4 | Accepted | 1,088 |
| 16 | 8 | Refused at existing action cap | 4,096 |

These probes are not larger-system energy certificates or equal-accuracy comparisons. The sixteen-site fourth-moment vector already has 8,256 orbit amplitudes. The full valence singlet basis also has Catalan growth. Ordinary molecular Hamiltonians generally fail the exact particle–hole gate and are not covered by this quotient.

The local Clifford frame and a connected-interval scalar moment compiler are now implemented; see [the connected dimer report](marginal_connected_dimer_moments.md). Exact eighth-order moments of a single product boundary use clusters of at most ten sites at arbitrary even chain length within the declared input bound. That new result does not replace the all-boundary, order-24 contraction above. Fixed-degree scalar Krylov energy bounds provably lose energy-density accuracy as chain size grows. General physical-marginal representability and the requested universal chemistry compiler remain unproved.

## Artifacts

- `experiments/marginal_symmetry_moments.py`: exact signed-orbit actions, moments, and renewal recurrence.
- `results/marginal_graded_hubbard8/discovery/singlet_moment_energy.py`: bounded energy replay using recomputed moments.
- `results/marginal_graded_hubbard8/singlet_moment_sharp/`: fresh range-two gap energy certificate and independent replay.
- `results/marginal_graded_hubbard8/singlet_moment_sharp_R3/`: slightly sharper range-three variant.
- `results/marginal_graded_hubbard8/symmetry_moments/independent_car_comparison.json`: exact H8 cross-check.
- `results/marginal_graded_hubbard8/symmetry_moment_scaling/receipt.json`: larger-chain successes and cap refusals.

```sh
python -S results/marginal_graded_hubbard8/discovery/singlet_moment_energy.py results/marginal_graded_hubbard8/singlet_moment_sharp/certificate.json
python -S -m unittest tests.test_marginal_symmetry_moments tests.test_marginal_singlet_moment_energy -v
```


Final validation: all **394 tests passed in 280.033 seconds**. The eight focused new tests passed in 10.835 seconds. Both sharp certificates passed separate `python -S` replay. The complete artifacts, including the retained embedding diagnostic data, are 99,080 bytes for the fresh range-two variant and 101,002 bytes for range three; the polynomial coefficient counts are not the entire certificate size.
