# Distilling the quartic proof into 36 square orbits

The accepted matched-model certificate can be distilled into **36 symmetry orbits of squared polynomials**, comprising 2,708 seed monomial coefficients. The rational orbit certificate is approximately 340 KB, compared with approximately 1.8 MB for the fuller factor certificate.

Its exact interval is

```text
3.280996277483316 <= E0 <= 3.281006695574010
width = 1.0418090694048992e-5
```

This is a size/accuracy tradeoff. The fuller Reynolds certificate has the tighter width 1.064960764376076e-6.

## Construction and acceptance

`experiments/marginal_distill.py` takes the accepted full SDP certificate, selects representatives of charge orbits, and extracts 562 candidate square directions from their integer factors. A linear program fits nonnegative weights plus the number-ideal multiplier. Ray normalization and a weight cap of 100 prevent the most severe numerical cancellation. The LP selected 36 rays and took about 1.18 seconds, excluding candidate construction and exact replay.

For each seed polynomial p and nonnegative rational weight a, the separate checker uses

\[
a\,\mathcal R(p^\dagger p)
=\frac{a}{|G|}\sum_{g\in G}g(p)^\dagger g(p)\succeq0.
\]

The checker validates fixed-charge seeds, rational nonnegative weights, orbital permutations, Hermitian number-conserving H and X, and the final residual. It does not require the supplied permutations to be Hamiltonian symmetries: positivity holds for each mapped square and the full coefficient residual is still checked. The final seed coefficients are rounded to a common denominator of 10^14; all resulting errors are included in the exact bound.

Two earlier LP proposals are retained in `results/marginal_distill` and its `normalized` subdirectory. Their numerical objectives exceeded the physical reference because large weights and cancellation produced actual coefficient mismatches of about 0.02. They are not accepted tight certificates. Ray normalization alone did not solve this issue. The bounded run's maximum numerical mismatch was 1.40e-7; exact rational replay, not LP status, determined the final 1.04e-5 interval.

## Remaining selection problem

This procedure compresses a certificate after a full quartic solve. It has not found the 36 useful directions directly from the Hamiltonian. The separate adaptive negative-eigenvalue selector failed to improve the energy in its seven-step test. The open algorithmic target is to discover a small family of energy-relevant positive squares without first solving the complete hierarchy.

## Replay

```sh
python3 -S -m experiments.marginal_orbit_certificate results/marginal_distill/compact_certificate.json
```

The recorded upper endpoint is recomputed from the integer symmetric amplitudes after exact comparison of H with the matched model. The full receipt is in `results/marginal_distill/compact_receipt.json`.

## Subsequent direct-discovery result

The older full-SDP compression above is retained. A separate pipeline now starts from cubic constraints, discovers scalar quartic directions through restricted dual pricing, and compresses the resulting certificate to 21 square orbits at certified width 1.6240708954489584e-6. No full quartic certificate supplies those directions. See `research/marginal_energy_selector_results.md`.
