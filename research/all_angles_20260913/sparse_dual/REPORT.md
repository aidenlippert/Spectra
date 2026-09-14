# Sparse dual / chordal probe

The bounded experiment uses the Hamiltonian-only active-space H4 fixture
(SHA-256 `8e64505deb4cc06a87f75e669b0f3073d15ebca72daa6d086e06eeea0b0f5120`),
not a prior certificate. It ranks two-word quadratic atoms by their exact H
overlap, solves a residual LP, exports integer factors, and replays them with
the exact CAR verifier. The sparse 64-atom run produced 21 blocks, 25 factor
nonzeros, 21 unique words, and 4 interaction edges. Its certified lower bound
is `-14.023370436342` with residual L1 `10.425057205342`, so it fails the
useful accuracy target. This is a real falsifying result for the small sparse
budget, not a copied certificate.

For comparison, the larger512-atomLP direct quadratic run produced 61
blocks, 122 factor nonzeros, lower bound `-8.408379706756`, and residual L1
`3.966127325924`. Both runs used the same H-only discovery and exact replay.
This is a larger finite atom family, not the full quadratic PSD/SOS cone.

The verifier rejected a deliberately malformed factor (`ValueError`). The
active-space paired upper witness is `-3.525719766889 Ha` from
`active_space_ladder/h4/upper.json` (SHA-256
`823b00624c2bad6ddd9be0b66ecbe6890435352934f895a458146a6596af4562`). These
large gaps show that this sparse family does not yet certify the target. All
pricing and omitted-family optimality costs remain explicit limitations.

Coverage: routes 1–8 (dual SOS, block pricing, coordinate/operator changes,
and adaptive support) and 17–20 (chordal/local decomposition). The probe
does not yet establish favorable asymptotic scaling. No complete asymptotic
treewidth test was performed. The next decisive test is
to run the same adapter against H6 and H8, recording complete candidate scans,
map nonzeros, wall time, and whether the exact verifier still accepts.

Scalable second-stage command (same H-only discovery, suitable for Lambda):

```sh
PYTHONPATH=/Users/aidenlippert/Documents/Spectra \
/opt/homebrew/Caskroom/miniconda/base/bin/python \
research/all_angles_20260913/sparse_dual/sparse_chordal_probe.py
```
