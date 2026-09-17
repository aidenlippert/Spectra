# Spectra: exact certificates for driven-model control

Read REPORT.md and DERIVATION.md before using the claim. This is a finite H8 control calculation from a specified rational MPS, not a certificate about an unknown exact ground state or a laboratory pulse.

## Complete independent acceptance

Python standard library only:

```sh
python3 -B -S replay.py --out ../fresh-control-replay
```

The output path must be new and outside this folder. Input files and old receipts are not modified. The replay reconstructs the actual molecular actions and the MPS, validates the reduced continuous-time paths, checks robust endpoint targets, and confirms the deliberate eight-coordinate refusal. It enumerates the 4,900-state balanced sector; that dependency is not concealed.

## Numerical proposal reproduction

Requires preinstalled NumPy, SciPy, and Numba; installs nothing:

```sh
python3 rediscover.py --out ../fresh-control-discovery
python3 -B -S ../fresh-control-discovery/code/exact_control.py --out ../fresh-control-discovery-replay
```

This repeats the H8-specific numerical procedure from the provided Hamiltonian and MPS, not from integrals. It includes failed/global-basis diagnostic phases and full-sector snapshot calculations. Numerical outputs are proposals until the exact second command accepts them. Software-version differences can change discovery; the exact accepting criterion remains the authority.

## Files

- inputs/: original Hamiltonian/MPS and rational basis/trajectory certificates.
- code/exact_control.py: full standard-library accepting implementation.
- code/test_exact_control.py: 20 algebra, independent-oracle, and refusal tests.
- code/*: numerical proposers and reference diagnostics.
- receipts/: accepted rational endpoints and timing evidence.
- development/: retained logs and numerical validation, including negative results.
- ACCOUNTING.json: separate computation stages and unavailable costs.
- RESULT.json: exact results and nonclaims.

The theorem is general finite-dimensional unitary/noisy-error mathematics. The molecular wrapper and proposed bases are specialized to the declared H8 test. Neither it nor the observed small reduced matrices demonstrate a generally inexpensive many-body constructor.
