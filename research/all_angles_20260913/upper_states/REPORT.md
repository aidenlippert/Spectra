# Upper-state route campaign

This bounded experiment covers routes 25–32 (variational states), 65–68 (stochastic methods), and 69–72 (embedding/coupled-cluster families). The implemented route is a Hamiltonian-only CISD witness: the lowest-diagonal determinant is selected from the frozen fixture H, all determinants up to two particle-hole substitutions are formed, the projected H is diagonalized, and the lowest eigenvector is rounded to an integer witness. No saved FCI coefficient is read during discovery. Each frozen witness is replayed against the original rational H by `streaming_reference_upper.py`.

| system | CISD basis | nonzero witness | upper (Ha) | lower from existing certificate (Ha) | interval width (Ha) |
|---|---:|---:|---:|---:|---:|
| H4 | 53 | 15 | -3.6573313621 | -3.66700010897 | 0.0096687469 |
| H6 | 262 | 60 | -6.3046674663 | -6.33311362063 | 0.0284461543 |
| H8 | 849 | 185 | -9.2025215477 | -9.25605501419 | 0.0535334665 |

The upper witness is valid and exact after rationalization, but these intervals do not meet the 0.0016 Ha target. The result is still useful: discovery takes about 0.02 s, 0.23 s, and 2.57 s locally (BLAS threads=1), while exact replay takes 0.004 s, 0.023 s, and 0.097 s. It demonstrates a cheap reproducible upper path whose cost is polynomial in the CISD subspace, with quality that degrades along this ladder.

Routes 65–68 (VMC, DMC, AFQMC, FCIQMC) remain prerequisite designs: their usual estimates require a rigorous bias/variance and fixed-node or phaseless error certificate before they can be upper endpoints. Routes 69–72 (embedding, DMFT, Green-function/self-energy, coupled-cluster) likewise need a certified truncation or residual bound. Routes 29–32 (PEPS/MERA/neural/circuit states) are plausible replacements for CISD, but no local implementation was claimed within the five-minute CPU bound. A natural next test is adaptive determinant selection by residual magnitude, then MPS/DMRG if a compatible dependency is available.

Failure checks: the replay API rejects non-integer amplitudes, duplicate/invalid states, empty witnesses, and witnesses outside the fixed-N sector. Changing a witness state to a wrong-sector bitstring or changing an amplitude to a float raises `ValueError` before energy production. All artifacts are in `results/all_angles_20260913/upper_states/`.

## Selected-CI refinement

The iterative H8 ladder (six residual rounds, 200 additions per round, capped at 2,000 determinants) reached `-9.2498074656 Ha` with 1,336 nonzero amplitudes. Discovery took 31.72 s and exact replay 0.48 s; the frozen projected matrix had 4,000,000 entries and 9,799 scored candidates. Against the existing lower, the width is about 0.00625 Ha, still failing the target. The full-sector legacy diagonal-reference scan touched 12,870 determinants and is recorded in the receipt.

`selected_ci.py` scores determinants coupled to the CISD seed using a residual/diagonal-denominator proxy, retains a bounded set, and re-diagonalizes. H6 reached `-6.3320509382 Ha` with a 393 determinant basis and 191 nonzero amplitudes in 0.50 s; H8 reached `-9.2258443581 Ha` with a 1,049 determinant basis and 385 nonzero amplitudes in 4.84 s. These improve the CISD upper by 0.0274 and 0.0233 Ha. The H6 gap to the existing lower is about 0.00106 Ha and H8 about 0.0302 Ha. Dense projected diagonalization is a bounded experiment and is not claimed scalable.
