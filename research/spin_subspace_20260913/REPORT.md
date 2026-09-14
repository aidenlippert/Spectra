# Spin coupling improves the H6 bound, but this operator family remains insufficient

**Thirty-two learned spin combinations reduce H6's certified interval from
13.010569 to 11.429079 mHa**, a 12.16% reduction. The same combinations used
separately give 12.108210 mHa, so keeping their cross terms improves the
found bound by another **0.679131 mHa**. Spin-sensitive interference supplies
useful missing information in this experiment.

The 1.6 mHa target is still unmet. The full spin-resolved pattern space gives
a 10.005183 mHa interval, and an exact dual proof establishes a **9.925834 mHa
physical error floor for this entire certificate family**. That proof covers
every subspace that could be selected from these generators. More directions
or more optimization inside the same family cannot reach the target.

The original compression question remains open. The ten-pattern Hamiltonian
and its collectively certified 0.108331264 mHa remainder are preserved. This
pass captures part of the hard joint correlations, but does not compress them
into a small sufficient description.

## The bounded experiment

The Hamiltonian is the same frozen rational electronic H6 chain: twelve spin
orbitals, six electrons, STO-3G, 1.4 Å spacing. Each retained spatial density
pattern Q_k is split into spin-up and spin-down components. Products
Q_k,spin a_i and a_i form four symmetry classes with 63 generators each.
Discovery uses neither a many-electron wavefunction nor old proof factors.

The common baseline is the complete quadratic fermionic cone. The previous
spin-summed joint density blocks are omitted because their exact ceiling is
essentially the baseline value; the ten-pattern Hamiltonian is unchanged.
The new learned combinations share a small joint anticommutator Gram block
within each symmetry class. All sixth-degree CAR terms cancel exactly;
none are discarded.

At each of four enrichment rounds, the current numerical dual proposes two
independent negative directions per symmetry class. Bundles adding four and
eight directions are both solved and exactly replayed. The continuation uses
certified gain per measured trial cost, including pricing. Both candidate
costs count toward the 180-second adaptive budget. Every round selected the
eight-direction bundle. Thus this run does **not** show that the selection
policy beats a fixed eight-direction schedule.

| Retained correlation combinations | Certified interval, mHa | Added joint Gram blocks |
|---|---:|---|
| 0: fresh quadratic baseline | 13.010569 | none |
| 8 | 12.765367 | 4 × 2 |
| 16 | 12.518908 | 4 × 4 |
| 24 | 11.910365 | 4 × 6 |
| **32** | **11.429079** | **4 × 8** |
| Same 32, separate contributions | 12.108210 | diagonal entries only |
| Full frame: 252 generators | 10.005183 | 4 × 63 |

The table shows the selected path and controls. Four additional candidate
solves at 4, 12, 20, and 28 combinations are preserved and charged in the
[campaign receipts](/Users/aidenlippert/Documents/Spectra/results/spin_subspace_20260913/campaign/adaptive/summary.json).
All eleven completed lower certificates, including the corrected separate
control, passed independent exact replay.

The 32-combination lower exceeds the **previous entire spin-summed cone's
exact ceiling** by 1.581045 mHa. Its improvement therefore cannot be explained
solely by optimizing that previous cone more accurately. The independent
control demonstrates a further measured gain from coupling the new
combinations; it is not an exact optimum claim for the separate cone.

## Compression and cost

The small subspace uses 32 combinations out of a 252-generator frame and
recovers **52.62% of the bound improvement achieved by the full-frame control**.
Its largest added Gram block has dimension eight; the unchanged quadratic
baseline still contains blocks up to dimension 36.

| Quantity | Learned 32-combination result | Full-frame control |
|---|---:|---:|
| Gram entries including baseline | 5,980 | 21,600 |
| Numerical coefficient-map nonzeros | 138,516 | 1,489,680 |
| Compact lower certificate | 67,798 B | 196,378 B |
| Expanded SOS certificate | 546,905 B | 1,312,728 B |
| Complete discovery/control time | **63.486 s** | **40.786 s** |

There is **no demonstrated end-to-end speed advantage**: discovering the
smaller subspace took longer than the full-frame control and produced a
weaker bound. The smaller proof has fewer Gram entries and fewer bytes, but
those savings do not remove candidate discovery and verification costs.

The 32 directions contain 2,016 nonzero coefficients over the existing spin
frames. Its certificate also contains 4,423 nonzero factor coefficients,
including the quadratic baseline. The ten spatial Hamiltonian patterns are
separate from these 32 learned correlation combinations. A small count of
either is not a count of physical configurations or a scaling theorem.

The final selected trial took 0.152 seconds for its additional construction,
1.896 seconds to solve, 0.028 seconds to export, and 5.651 seconds for exact
acceptance. The 63.486-second total also includes shared initial construction,
the baseline, pricing, all eight candidate trials, imports, and writes.
The search stopped at its four-round limit before exhausting the time budget.

The first separate-control attempt failed after 5.235 seconds, before solving,
because its reporting code compared a relative path with an absolute root.
That log and source snapshot are preserved. The corrected control completed
in 11.997 seconds under a fresh 90-second watchdog. Its descriptors reuse the
63.486-second adaptive discovery; that ancestry remains charged. The initial
campaign took 109.511 seconds including the failed control and full-frame
run. The retry adds 11.997 seconds. These are single-run observations.

## An exact limit on the full spin frame

The full-frame numerical result was `optimal_inaccurate`. It was followed by
a separate, bounded exact diagnostic. Rational affine projection and a 10^-7
mixture with the uniform fixed-N trace functional produced an accepted dual
witness. Two smaller tested mixtures failed exact PSD checks.

The witness satisfies normalization, coefficient bounds, Hermiticity, all
79 one-body number-ideal identities, every quadratic moment PSD condition,
and all four complete 63-dimensional spin-frame PSD conditions. Its ceiling
on the molecular lower bound is

\[
C=-6.343039454711732\ldots\ \mathrm{Ha}.
\]

The [fresh obstruction audit](/Users/aidenlippert/Documents/Spectra/results/spin_subspace_20260913/audit.json)
also replays the existing full-cubic physical lower and its separate spectral
residual proof against exactly the same Hamiltonian and particle sector.
That lower is L*=-6.333113620625836... Ha. Consequently, every lower L from
the present full frame or any selected subspace obeys

\[
E_0(H6)-L\ge L_*-C
=9.925834085896\ldots\ \mathrm{mHa}>1.6\ \mathrm{mHa}.
\]

With the frozen reference upper, the narrowest possible interval in this
family is at least **9.980828 mHa**. The full-frame export is within
0.024355 mHa of that ceiling. The statement is about the specified operators,
one-body number ideal, coefficient-L1 residual rule, and fixed tail transfer.
It does not rule out compact representations using other joint operators.

## What should change next

The next pass should use this exact counterexample to identify necessary
constraints **outside the current density-times-annihilator frame**, while
keeping the ten-pattern Hamiltonian and its tail fixed. Spin-changing density
components and commutators between patterns are concrete candidates; neither
has been shown here to provide the required gain.

The useful question is now which few additional joint operators remove the
remaining nonphysical moment assignments at an affordable total cost. Further
selection inside the present 252-generator frame cannot answer the accuracy
target, regardless of how well it is optimized. The full-frame counterexample
provides a precise test for the next representation rather than an invitation
to keep tuning this one.

## Verification and reproducibility

- **Eight focused tests passed**, across the initial six-test run and two
  additional full-frame gate tests. They cover two-stage coordinate
  composition, exact adjoint pairing, tail transfer, invalid input refusal,
  numerical pricing versus rational Gram matrices, independence of proposed
  operators, complete reconstruction of the full frame, and rejection of the
  old spin-summed dual by the stronger spin gate.
- Eleven molecular lower certificates and the old exact cone comparison
  passed fresh `python -S` replay in **59.085 seconds**, including the
  existing rational upper witness. No NumPy, SciPy, CVXPY, or PySCF was loaded
  on that accepting path.
- The full-frame dual construction, exact replay, and export took
  **59.698 seconds** under a separate 120-second watchdog. Its fresh combined
  audit took **37.624 seconds**, including 31.538 seconds for the new dual
  and 6.035 seconds for the independent physical interval proof.
- All **110 files** in the preceding pass's manifest still match their hashes.
  Its implementation, fixtures, patterns, certificates, and reports are
  unchanged. This pass writes only its new research and result directories.

The old rational upper uses 200 nonzero determinant amplitudes. Its prior
exponential wavefunction discovery is not part of the compact lower-bound
method. No full physical sector is enumerated in the new discovery or
acceptance. These are finite-basis electronic-Hamiltonian statements; they do
not establish general many-body solvability or experimental accuracy.

The [proof](/Users/aidenlippert/Documents/Spectra/research/spin_subspace_20260913/PROOF.md)
describes the certificate and dual implications. All logs, rejected proposals,
costs, and exact witnesses are in
[results/spin_subspace_20260913](/Users/aidenlippert/Documents/Spectra/results/spin_subspace_20260913).

From the repository root, using the configured Python environment:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest \
  research.spin_subspace_20260913.test_core \
  research.spin_subspace_20260913.test_discovery -v
python -S -m research.spin_subspace_20260913.replay \
  --campaign results/spin_subspace_20260913/campaign \
  --separate results/spin_subspace_20260913/separate_retry \
  --out /tmp/spin-subspace-replay.json
python -S -m research.spin_subspace_20260913.audit \
  --results results/spin_subspace_20260913 \
  --out /tmp/spin-subspace-audit.json
```

`campaign.py --out <new-directory>` repeats the corrected adaptive workflow and
both controls, preserving the budgets. `full_dual.py` implements the separate
post-campaign diagnostic. Existing output directories are refused by discovery.
