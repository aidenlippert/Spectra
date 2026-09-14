# H6: ten collective interaction patterns with a certified remainder

The frozen H6 Hamiltonian now has a representation using **ten collective
density-square interaction patterns**, an explicit one-body operator, and a
remainder whose entire operator interval is **0.108331264 mHa** wide. This is
below the 1.6 mHa target used to test the remainder. The patterns are extracted
from H6's two-electron coefficients; no ground-state wavefunction or previously
discovered proof factors are used in their discovery.

The representation handles the molecular remainder that the previous
four-orbital/collective-charge model left unresolved. Its ten patterns act
across all six spatial orbitals. Solving their joint many-electron behavior
remains necessary: this pass proves useful compression of the Hamiltonian
and its remainder, not a tight new ground-energy interval from a compact
retained-system solver.

## Representation and the collective rule

Each pattern has the form A_k=Σ_pqσ L^k_pq a†_pσ a_qσ, with real symmetric
coefficients. The retained model contains ½Σ_k λ_k A_k² plus its one-body and
constant terms. The square retains interference between all its components.
Exact CAR expansion reconstructs the contractions and charges the remaining
rounding mismatch explicitly. For H6 that mismatch is bounded by
3.2×10⁻¹¹ Ha.

Before selecting factors, the conserved total-number direction is removed
exactly and its couplings are folded into the one-body operator. The discarded
coefficient matrix is then bounded by two multiples of a fixed feature
metric. A particle-number/spin Casimir identity bounds all those features
together, without visiting physical configurations. Both coefficient-matrix
PSD tests use exact rational arithmetic.

For the exported ten-pattern model H₁₀, the precise H6 operator statement is

\[
H_{10}-3.632\times10^{-9}\,I
\preceq H_{\mathrm{H6}}
\preceq H_{10}+1.08327632\times10^{-4}\,I
\quad\text{(Hartree, fixed }N=6\text{)}.
\]

This holds for every state in the sector. Any certified interval for H₁₀
can therefore be transferred to H6 with only this additional width. The
[proof](/Users/aidenlippert/Documents/Spectra/research/molecular_collective_20260913/PROOF.md)
derives the number identity, collective bound, factor acceptance, and the
rank and energy obstructions.

## Size and cost

The inputs are straight hydrogen chains at 1.4 Å spacing in the frozen
STO-3G/RHF orbital fixtures. Energies are electronic Hartree; these bounds do
not include basis-set or experimental-model error.

| Fixture | Spin orbitals / electrons | Retained patterns | Nonzero pattern coefficients | Certified remainder width | Tail certificate | Exact tail replay |
|---|---:|---:|---:|---:|---:|---:|
| H4 | 8 / 4 | 6 | 30 | 0.048883 mHa | 1,352 B | 0.0154 s |
| H6 | 12 / 6 | 10 | 105 | 0.108331 mHa | 3,667 B | 0.0753 s |
| H8 | 16 / 8 | 14 | 252 | 0.187166 mHa | 7,903 B | 0.3156 s |
| H10 | 20 / 10 | 18 | 495 | 0.286368 mHa | 14,827 B | 0.8193 s |

These are minimum ranks for a 1.6 mHa remainder budget **within this particular
coefficient-matrix envelope and collective scalar bound**. Exact positive
subspace witnesses prove that one fewer factor cannot meet the budget, for
any orientation of the retained rank-limited matrix. Passing constructions
match these lower bounds. This does not exclude another operator bound or a
different pattern representation.

The observed counts are 2s−2 for s=4,6,8,10 spatial orbitals. Four cases are
not an asymptotic scaling theorem. Each factor still needs up to s(s+1)/2
coefficients, and the two PSD checks have that dimension: 10, 21, 36, and 55.
Reading and reconstructing the original Hamiltonian remains part of the cost.

Number elimination matters. With ten factors but without eliminating the
number direction, the H6 remainder bound is 127.387790 mHa. Eliminating that
direction brings it to 0.108331 mHa. Nine centered factors still give
104.250650 mHa using the implemented bound; the exact rank witness excludes
reaching 1.6 mHa with rank nine in this grammar.

The H6 retained model, including its one-body operator and constant, occupies
6,683 bytes in compact JSON, or 13,560 bytes as saved indented JSON. Its tail
certificate is separate and requires the original Hamiltonian to replay.
The source H6 input has 918 canonical terms. The model is available as
[retained_model.json](/Users/aidenlippert/Documents/Spectra/results/molecular_collective_20260913/campaign/h6/retained_model.json).

The complete measured campaign took **23.528 seconds**. It includes 45
centered rank trials, four uncentered controls, four rank obstructions, eight
supporting-bound energy calculations, exact replays, and artifact writes.
H6's complete portion took 3.681 seconds. Its initial extraction took 0.0971
seconds; the selected factor proposal took 0.0190 seconds, followed by the
0.0753-second accepting replay. These are single-run measurements.

An initial profile found repeated copying during CAR accumulation dominated
larger-input extraction. Accumulating coefficients directly reduced the H10
extraction from about eleven seconds in that exploratory run to 0.764 seconds
in the recorded campaign. The algebra and acceptance conditions are unchanged.
Exploratory runs, interpreter startup, earlier fixture/integral generation,
and development time are outside the campaign timer.

## What still makes the ten-factor model difficult

We tested the inexpensive lower-bound family obtained from
A_k²≥2c_k A_k−c_k². For any choice of centers c_k, the resulting one-body
Hamiltonian is treated jointly, with an exact PSD certificate for its
fixed-number lower bound. A convex program proposes the centers. This produces
valid energy intervals using the HF determinant as the upper witness.

For the ten-factor H6 model, the resulting lower bound on the original H6
Hamiltonian is −6.86518608479 Ha. The HF upper is −6.125774106371 Ha, giving a
wide interval. The physically validated reference ground energy is about
−6.33305862623 Ha.

The weakness is not an unfinished optimization. An exact rational density
matrix with 0≤γ≤I and tr γ=6 proves that **every choice of the square supports**
has a lower-bound ceiling of −6.86518596575 Ha. Our found lower is within
1.191×10⁻⁷ Ha of that ceiling. An independently replayed earlier lower
certificate proves E₀(H6)≥−6.345083962727 Ha. Together they establish

\[
E_0(H6)-\text{any lower bound in this supporting family}
\ge 0.520102003019\ \mathrm{Ha}.
\]

Keeping all available nonzero rounded factors still leaves a certified gap
of at least 0.520100942364 Ha for the same supporting family. More factors or
more center optimization therefore cannot repair this particular solver.
The exact obstruction is replayed by `audit.py`, with hashes of both its new
dual witness and the independent prior lower proof.

The retained patterns also have explicitly nonzero commutators. Treating all
their fluctuations as independently removable loses a substantial part of
the fermionic constraints. The obstruction concerns the separate supporting
planes; it says nothing against a solver using joint constraints on the
patterns.

## Validation and receipts

- **14 exact tests passed.** They cover the general CAR ingredients on small
  systems, number elimination on several particle sectors, matrix bounds,
  rank witnesses, dual feasibility, and rejection of corrupt or out-of-budget
  inputs.
- **49 tail certificates, eight supporting-bound certificates, and four rank
  obstructions** passed a fresh `python -S` process in 14.575 seconds. It
  imported no NumPy, SciPy, CVXPY, or PySCF. An earlier attempt failed in the
  reporting wrapper because a relative path was compared with an absolute
  root; its log is preserved. The corrected complete replay passed.
- Full physical-sector controls use 70 H4 and 924 H6 states, separately from
  discovery and acceptance. Their complete difference-operator spectra lie
  within the certified envelopes. Numerical ground-energy changes are
  **0.002290 mHa for H4** and **0.004713 mHa for H6**. These smaller observed
  errors do not replace the rigorous worst-state bounds.
- The full controls took 2.089 seconds and 2,585,324 word/state actions,
  including exact replays of existing rational reference upper witnesses.
  Their exponentially growing physical-sector work is disclosed and capped
  at twelve spin orbitals.
- Exporting the retained models and replaying the exact energy obstruction
  took 10.860 seconds, including 5.121 seconds to verify the independent prior
  lower proof. That proof's earlier discovery cost is not part of the new
  factor-discovery performance claim.

All data, certificates, model exports, source hashes, and logs are in
[results/molecular_collective_20260913](/Users/aidenlippert/Documents/Spectra/results/molecular_collective_20260913).
Existing solvers, prior certificates, and fixtures were preserved.

## Next course

Keep the ten H6 patterns and their certified tail. The next bounded test
should impose **joint fermionic constraints on those patterns**, with all
cross terms retained. Products of the learned density operators with fermion
creation/annihilation operators give Hamiltonian-derived candidates for such
constraints; their CAR identities and construction cost must be checked.

The success criterion is a useful lower bound for the retained model plus
the already controlled remainder, with pattern description, construction,
solve, and replay costs all charged. This directly continues the original
compression question. The exact supporting-plane obstruction makes another
round of independent square optimization unnecessary.

Density-square and nested low-rank Hamiltonian factorizations have established
antecedents, including [Motta et al., *Low rank representations for quantum
simulation of electronic structure*](https://arxiv.org/abs/1808.02625). That work
targets quantum simulation primitives; it does not establish the ground-energy
compression claimed or tested here. This pass supplies exact, scoped
certificates around this representation in the present project. No first-in-
literature or general many-body-solvability claim is made.

From the repository root, reproduce the checks with:

```sh
python -S -m unittest research.molecular_collective_20260913.test_core -v
python -S -m research.molecular_collective_20260913.replay \
  --campaign results/molecular_collective_20260913/campaign \
  --out /tmp/molecular-collective-replay.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  -m research.molecular_collective_20260913.validation \
  --campaign results/molecular_collective_20260913/campaign \
  --out /tmp/molecular-collective-validation.json
python -S -m research.molecular_collective_20260913.audit \
  --campaign results/molecular_collective_20260913/campaign \
  --out /tmp/molecular-collective-audit.json
```

`campaign.py --out <new-directory>` reruns discovery and refuses an existing
campaign directory. The final model export step intentionally rewrites only
the exports created by this pass.
