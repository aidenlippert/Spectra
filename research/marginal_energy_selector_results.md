# Direct discovery of a compact quartic certificate

The matched ten-mode, five-particle test now has a **directly discovered** certificate with interval width **8.157088389144001e-7**. A subsequent linear-programming compression reduces it to **21 square orbits**, with width **1.6240708954489584e-6**. Neither step takes directions from a full quartic SDP solution or a many-body eigenvector.

The model uses hopping 1/5 and density interaction coefficient one. All widths below are in those model energy units.

## Discovery algorithm

The restricted optimization begins with the complete cubic SOS representation, reduced to 64 PSD blocks (largest dimension 6; 311 PSD scalar variables). It uses the 39 invariant coefficient equations and 19 invariant number-multiplier directions described in the dimension note.

For the current restricted solution, each of the 30 representative quartic dictionaries is priced against the equality dual. A negative eigenvector of the symmetrized adjoint Gram map proposes an orbit of a positive square. Only a nonnegative scalar weight for that direction is added to the optimization. Each scalar weight is capped at 100, an explicit extra restriction on the search. This cap can weaken discovery; it does not weaken certificate validity.

Three selectors were compared with the same 117-cut budget:

| Selector | Rounds | Added scalar square directions | Exact interval width |
|---|---:|---:|---:|
| One-step energy lookahead, with batch fallback | 51 | 117 | 9.928164433493868e-7 |
| Most negative direction, one at a time | 117 | 117 | **8.157088389144001e-7** |
| Most negative directions in batches of up to eight | 15 | 117 | 1.0306578739844184e-3 |

The lookahead run took 26.45 seconds in total, including its assembly, searches, and certificate checks. The single-direction run took 27.85 seconds and the batch run 12.14 seconds. These are individual observed timings, not controlled repeated performance measurements; the two baseline processes ran concurrently. The comparison establishes **no advantage for energy lookahead on this instance**. Simple spectral pricing also works after sufficient iterations in this symmetry-reduced formulation.

This corrects the interpretation of the earlier seven-step failure. That bounded earlier run did not improve the bound, but it did not prove that the selection principle is incapable of doing so. The new experiments also change the symmetry formulation and iteration budget, so their success does not isolate a single causal improvement.

The lookahead run stopped because no eigenvalue below -1e-7 was detected in the priced dictionaries. This is a numerical stopping condition, not a rational proof of full moment feasibility or exact convergence of the hierarchy. The single-direction baseline stopped at its cut budget. Their independently checked lower/upper intervals are the rigorous outcomes.

Column generation for PSD and SOS approximations is established methodology; see [Ahmadi, Dash, and Hall](https://arxiv.org/abs/1512.05402). The result here is this fermionic implementation, comparison, and exact certificate—not a claim to have invented column generation.

## Exact certificate and scalar normalization

The single-direction certificate gives

```text
3.2810058798651713 <= E0 <= 3.281006695574010
```

The numerical search maximizes b minus a coefficient residual norm. Its raw export contained an approximately -25.387685 scalar residual, offset by a correspondingly large b. This is largely an arbitrary scalar offset, not a 25-unit error in the energy bound. Absorbing that scalar into b **exactly** leaves the certified lower endpoint unchanged. The normalized residual L1 is **7.576743268515014e-8**. The normalization is separately tested and replayed; original artifacts are retained.

At the final step the optimization contains 311 cubic PSD scalar variables plus 117 added scalar weights, compared with 2,082 PSD scalar variables in the full symmetry-reduced quartic SDP. This is a reduction in optimized cone variables. Candidate construction still enumerates complete quartic word dictionaries, and pricing still diagonalizes candidate moment matrices up to dimension 186. No general cost scaling improvement follows from the variable count alone.

A subsequent optional `--pricing stabilizer` path replaces those full eigensolves by batched subblock eigensolves of maximum dimension 14. It preserves the priced spectrum in the checked instances and reduces isolated pricing time by about fivefold. Full-map preparation remains, and the adaptive cut trajectory changes with degenerate eigenvector choices. See [the separate benchmark and complete-run result](marginal_reduced_pricing_results.md); the table above remains the original full-pricing comparison.

## Compression after direct discovery

The resulting 258 square seeds, including factors from the cubic base, were passed to a linear program. It selected 21 orbit seeds containing 2,044 monomial coefficients. The LP took approximately 0.08 seconds, excluding candidate-map construction and exact replay. Common-denominator rational rounding and exact replay give

```text
3.281005071503115 <= E0 <= 3.281006695574010
width = 1.6240708954489584e-6
```

The complete provenance chain is therefore **cubic base → restricted moment pricing → scalar quartic additions → LP compression → rational verification**. The older full-quartic certificate is not an input to this chain.

## Reproduction and replay

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_energy_selector --strategy negative --iterations 117 --max-cuts 117
python3 -S -m experiments.marginal_orbit_certificate results/marginal_energy_selector/m10_steps117_negative_cap117/normalized_certificate.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_distill --source results/marginal_energy_selector/m10_steps117_negative_cap117/normalized_certificate.json --output-name marginal_distill_direct
python3 -S -m experiments.marginal_orbit_certificate results/marginal_distill_direct/compact_certificate.json
```

The upper endpoints are independently recomputed from integer symmetric amplitudes after exact Hamiltonian comparison. The full receipts include those endpoints. Generalization beyond the matched symmetry class, large-system pricing cost, required operator degree, and unrestricted chemical accuracy remain unproved.
