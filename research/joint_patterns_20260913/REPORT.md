# The ten-pattern Hamiltonian survives; its first joint constraint family fails

The H6 Hamiltonian still has **ten learned density-square patterns** and a
collectively certified remainder only **0.108331264 mHa** wide. This pass tested
whether a small set of joint fermionic constraints could solve that retained
model accurately. The answer for the tested family is **no**, with an exact
obstruction rather than an inconclusive numerical result.

The best new H6 interval is **13.010458 mHa**, against a 1.6 mHa target. An
exact dual certificate proves this family's intervals cannot become narrower
than **13.010125 mHa** with the fixed reference upper. A separately replayed
prior physical lower proves that every lower bound from this family misses
the true H6 ground energy by at least **12.955130 mHa**.

There is a concrete next candidate: resolving the same ten patterns into their
spin-up and spin-down components reveals an exact fermionic constraint that
the current family misses. This diagnoses missing joint information inside
the retained model. It does not establish a new tight energy bound.

## What was held fixed and what was tested

The inputs are the existing rational STO-3G hydrogen chains at 1.4 Å spacing:
H4 has eight spin orbitals and four electrons; H6 has twelve and six. These
are bounds on their electronic finite-basis Hamiltonians.

Every H6 trial uses the same ten-pattern Hamiltonian and the same tail
certificate. Only the constraints on that model change. For each learned
density Q_k, form Q_k a_i and add a_i. A joint block permits arbitrary linear
combinations B and uses {B^dagger,B}>=0, retaining every cross-pattern term.
The degree-six terms cancel by exact CAR identities, leaving coefficient
rows of degree at most four. A complete quadratic fermionic baseline and
a one-body multiplier of the particle-number identity remain in every case.

The separate-pattern ablation removes the cross-pattern Gram entries. Two-
and four-pattern trials restrict the extra constraints, not the Hamiltonian.
No ground-state vector or prior proof factor is a discovery input.

| Case | Certified interval, mHa | End-to-end time | Added Gram blocks | Compact lower proof |
|---|---:|---:|---|---:|
| H4 quadratic baseline | 4.482281 | 1.323 s | none | 11,466 B |
| H4 joint six patterns | 4.478848 | 1.725 s | 4 × 14 | 15,460 B |
| H6 quadratic baseline | 13.010569 | 2.744 s | none | 40,633 B |
| H6 separate ten patterns | 13.015592 | 27.777 s | 40 × 6 | 52,041 B |
| H6 joint two patterns | **13.010458** | 11.606 s | 4 × 9 | 42,346 B |
| H6 joint four patterns | 13.014253 | 11.832 s | 4 × 15 | 48,612 B |
| H6 joint ten patterns | 13.012080 | 54.013 s | 4 × 33 | 65,072 B |

All seven cases completed inside the predetermined 90-second budget, including
process startup, construction, solving, rational export, exact acceptance,
and file writes. The campaign took **111.030 seconds**. The slightly worse
four- and ten-pattern exports do not reverse cone inclusion: their PSD-clipping
and residual penalties differ. The accepted dual supplies the stronger answer
about the whole ten-pattern cone.

This pass does not set a new best H6 energy bound. The earlier full cubic
precision proof already achieves a **0.054994 mHa** interval with a much larger
discovery dictionary and separate residual factors. The present experiment
tests whether the compressed Hamiltonian admits a much smaller sufficient
constraint set.

## Why more optimization of this family will not help

The exact dual functional obeys all 79 one-body number-ideal identities,
normalization, coefficient bounds, all 26 quadratic PSD conditions, and all
four joint 33-dimensional PSD conditions. Its ceiling on the transferred
molecular lower bound is **−6.346068751063810... Ha**. The best new lower is
within **0.000333 mHa** of that ceiling.

The [final audit](/Users/aidenlippert/Documents/Spectra/results/joint_patterns_20260913/audit_final.json)
also replays the old full cubic precision certificate and its spectral
residual proof against this exact molecular fixture. Its independent lower
is −6.333113620625836... Ha. Comparing these two exact certificates proves
the **12.955130 mHa physical error floor** for the tested family. This argument
does not rely on a numerical ground energy or on a solver status.

The bound applies to this particular anticommutator dictionary, number
multiplier, and coefficient-L1 residual rule. It does not refute compression
of many-electron correlations more generally. The [proof](/Users/aidenlippert/Documents/Spectra/research/joint_patterns_20260913/PROOF.md) gives
the certificate identity and the dual inequality.

## The missing information exposed by spin resolution

Splitting each Q_k into spin-up and spin-down components reuses all ten
spatial coefficient matrices. Four diagnostic blocks then have dimension 63.
The accepted dual violates a positive anticommutator condition from these
generators: its exactly checked normalized expectation is
**−0.011516517075511...**. This is a moment diagnostic, not a Hartree gain.

The exported [direction](/Users/aidenlippert/Documents/Spectra/results/joint_patterns_20260913/spin_diagnostic/separator.json)
has 63 nonzero coefficients and expands to 156 monomials. It gives a concrete
necessary condition that can rule out the current dual witness. Finding it,
including replaying the dual first and verifying the direction exactly, took
**12.997 seconds** under a separate 60-second budget.

The next bounded energy test should add a small set of these spin-resolved
combinations, keeping the ten spatial patterns and the same collective tail.
Selection should be justified by actual certified bound gain and complete
construction/solve/replay cost. A negative moment alone does not establish
that the resulting lower bound will improve enough; another feasible dual
could replace the one excluded here.

## Complete cost and verification

The joint ten-pattern H6 run took 4.625 seconds to construct its maps,
42.760 seconds to solve, 0.045 seconds to export, and 5.443 seconds for exact
acceptance; process startup and writes bring the total to 54.013 seconds.
There are 919 coefficient rows, 10,080 Gram entries including the quadratic
baseline, and 1,726,272 numerical map nonzeros. Construction visits 103,068
monomial word pairs and accounts for 113,595,048 dense outer-product entries.
Small Gram blocks do not make this coefficient work disappear.

The compact H6 proof contains 7,384 nonzero factor coefficients and expands
to 743,563 bytes in the existing SOS format. Its 65,072-byte compact size
excludes the frozen Hamiltonian and tail inputs needed to replay it. The
preceding retained-model descriptor is 6,683 bytes in compact JSON and its
tail certificate is 3,667 bytes; those remain separate artifacts.

- **16 focused tests pass.** They check independent fermion actions, cubic
  cancellation, cross-term maps and basis changes, tied factor rounding,
  rational Gram contraction, affine dual repair, spin decomposition, and
  rejection of corrupted bindings, invalid factors, false ideal/symmetry
  claims, and nonseparating directions.
  The 14 exact-only tests also pass under `python -S`.
- Seven lower certificates and both prior rational upper witnesses passed
  fresh `python -S` replay in **27.647 seconds**. No NumPy, SciPy, CVXPY, or
  PySCF was loaded on that accepting path.
- The exact dual was constructed and replayed in **14.269 seconds**, inside
  its separate 120-second budget. All 30 PSD checks passed. Its 10^-7 trace
  mixture and the two smaller rejected attempts are recorded.
- The final combined dual, prior physical proof, upper, and spin-direction
  audit passed in **13.452 seconds**, including **5.784 seconds** to replay
  the prior physical lower. That prior proof occupies 970,728 bytes plus
  95,108 bytes for its spectral residual certificate.
- All **155 files** in the preceding molecular pass's manifest still match
  their recorded hashes. Existing solvers, fixtures, reports, and certificates
  were preserved.

These are single-run timings, not a controlled speedup claim. No full physical
sector was enumerated in the new discovery or accepting paths. Existing
reference uppers use 20 H4 and 200 H6 nonzero determinant amplitudes; their
prior exponential wavefunction discovery is not removed by this method.
Development, interpreter/test runs outside the measured campaign, the extra
fresh dual replay, and an initial audit using a weaker old comparison proof
are outside the campaign timer and remain in the logs. The initial comparison
audit took 12.663 seconds; the final audit uses the stronger spectral proof.

The underlying partial three-particle positivity ideas have established
antecedents, including [Mazziotti, Phys. Rev. A 71, 062503 (2005)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.71.062503)
and [Phys. Rev. Lett. 117, 153001 (2016)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.117.153001).
This pass supplies a scoped test, exact obstruction, and explicit missing
constraint for the project's learned patterns; it makes no priority claim
or general many-body solvability claim.

Reproduce from the repository root with the configured Python environment:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest \
  research.joint_patterns_20260913.test_core \
  research.joint_patterns_20260913.test_discovery \
  research.joint_patterns_20260913.test_dual \
  research.joint_patterns_20260913.test_spin_diagnostic -v
python -S -m research.joint_patterns_20260913.replay \
  --campaign results/joint_patterns_20260913/campaign \
  --out /tmp/joint-pattern-replay.json
python -S -m research.joint_patterns_20260913.audit \
  --results results/joint_patterns_20260913 \
  --out /tmp/joint-pattern-audit.json
```

`campaign.py --out <new-directory>` repeats discovery under the outer watchdog
and refuses an existing campaign directory. `dual.py` and
`spin_diagnostic.py` expose the additional proposal and exact-replay entry
points. All receipts and logs are in
[results/joint_patterns_20260913](/Users/aidenlippert/Documents/Spectra/results/joint_patterns_20260913).
