# H4 molecular certificate

The fixture is rectangular H4 with coordinates (0,0,0), (1,0,0), (0,1.5,0), (1,1.5,0) Angstrom, in STO-3G. PySCF RHF canonical MOs produce 8 spin orbitals and the 4-electron FCI total energy is -2.1249032164909356 Ha (electronic energy -4.475896658361423 Ha; nuclear repulsion 2.350993441870487 Ha).

An independent 70-dimensional CAR diagonalization of the exported rational Hamiltonian gives -4.475896658361418 Ha, differing by 4.4e-15 Ha. The finite coefficient export has an exact coefficient l1 perturbation budget 6.0682e-14; basis, geometry, and integral algorithm errors are separate and are not certified.

The generic alpha/beta charge-separated degree-3 SOS search used 25 blocks, largest dimension 92, and ran 106.8 seconds with one thread. Exact replay gives lower -4.475905060890658 Ha and an independently rounded Rayleigh upper witness, width 3.19546e-5 Ha. The interval is for the rational CAR fixture; the upper witness is stored in `results/marginal_molecule/h4_degree3_certificate.json` and replayed by `experiments/marginal_transfer_verify.py`.

## Compact coefficient and upper-witness refinement

The original artifacts above are preserved. `experiments/marginal_molecule_accept.py` rounds the canonical coefficients to a common denominator of 10^12, reuses the original SOS factors, and regenerates an integer Rayleigh witness. Exact replay of `results/marginal_molecule/h4_compact_certificate.json` gives lower -4.47590506089335 Ha, upper -4.475896658362239 Ha, and width **8.40253111066967e-6 Ha**. The conservative coefficient perturbation enclosure, including the original rationalization, is 2.1582e-11 Ha. This encloses coefficient export error relative to computed floating integrals; it does not certify the underlying integral algorithm or finite-basis approximation.

Replay without numerical dependencies:

```sh
python3 -S -m experiments.marginal_transfer_verify results/marginal_molecule/h4_compact_certificate.json
```
