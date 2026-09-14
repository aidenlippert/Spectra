# Collective elimination around retained interference

This pass implements the original question in a restricted model family:
retain a small interacting system and its collective coupling patterns, while
certifying all other modes together. The construction works, including for
dispersed bath energies. It does **not** establish that the difficult part of
a molecular Hamiltonian naturally has this form.

The key distinction is that the remainder is now actually eliminated by a
short occupation rule. A large determinant list, coefficient map, or bath
matrix is not hidden behind the small accepting matrices. The family was
constructed to have the required collective structure; discovering suitable
structure in molecular two-electron interactions remains the central problem.

## What was built

Four active fermion modes have nontrivial hopping and two-body interactions.
They couple to positive and negative bath bands through two specified linear
combinations. Collective charging interactions are reduced exactly using
fixed total particle number. Each band has B = L² modes, giving M = 4+2B
physical modes and N = B+2 particles. Per-orbital hybridization scales as 1/L;
the total coupling strength is held fixed as the bath grows.

Bath particle and hole excitation energies are replaced by bin minima and
maxima, giving two operator envelopes with the same exact filled-bath energy,
active interactions, and coupling. Inside each envelope, one bright mode per
bin couples to the active modes. All dark occupations are minimized by their
energy and multiplicity while maintaining total particle number. A rational
CAR metric represents each bright mode without square roots in acceptance.

The certificate reconstructs every feasible retained particle sector, proves
the lower endpoint by exact PSD checks, and obtains the upper endpoint from
an exact rational Rayleigh witness or a proved excitation bound. Floating
eigenvectors are proposals only. The theorem and failure boundary are in
[PROOF.md](/Users/aidenlippert/Documents/Spectra/research/collective_interference_20260913/PROOF.md).

Here “pattern” means a retained collective bath combination, together with
its specified active coupling. This is a different representation from the
previous molecular local-cubic moment directions. It does not yet compress
the four-mode active interaction itself into fewer patterns.

## Size, accuracy, and actual cost

All energies in this table are **abstract model units**, not Hartree.

| Case | Physical modes / particles | Retained modes / bath patterns | Certified interval width | Largest proof matrix | Accepting replay |
|---|---:|---:|---:|---:|---:|
| Degenerate bands, L=64 | 8,196 / 4,098 | 6 / 2 | 8.164×10⁻¹² | 20×20 | 0.0212 s |
| Dispersed bands, L=64 | 8,196 / 4,098 | 6 / 2 | 0.01508851 | 20×20 | 0.0206 s |
| Dispersed, two bins per band, L=2 | 12 / 6 | 8 / 4 | 0.00503432 | 70×70 | 1.9144 s |
| Dispersed, two bins per band, L=16 | 516 / 258 | 8 / 4 | 0.00751474 | 70×70 | 2.1037 s |

For the one-bin dispersed family the interval width stayed 0.01508851 as M
grew from 12 to 8,196. This has a structural explanation: after subtracting
the exactly known filled-bath energy, the two normalized retained envelope
Hamiltonians are independent of B for B≥4. Both coupling and spectator
minimization remain bounded. The result is stronger than merely seeing a flat
timing curve on the tested sizes.

At L=16, refining each band from one bin to two reduced the interval by about
50.2%, while accepting replay increased from about 0.021 to 2.104 seconds.
The retained calculation remains exponential in its mode count. With one/two
bins it constructs at most 64/256 retained occupation labels and 924/12,870
matrix entries per envelope. It enforces the eight-retained-mode budget.
The 8,196-mode calculation eliminates 8,190 dark modes.

The largest one-bin certificate occupies 804 bytes for the dispersed case;
its model occupies 3,185 bytes as the saved indented JSON file, or 1,034 bytes
as compact JSON. The two-bin L=16 certificate occupies 2,020 bytes. These
sizes exclude the reusable checker and mathematical proof. Integer precision
and model/witness bit lengths grow with bath multiplicity; constant matrix
dimension is not a claim of literally constant bit cost.

The 22-case campaign took 17.467 seconds inside the campaign process: 9.043
seconds for construction/proposal/exact lower checking and 8.291 seconds for
accepting replay, with the rest in import initialization, orchestration, and
file output. A separate fresh `python -S` process replayed all 22 certificates
in 8.598 seconds without NumPy, SciPy, CVXPY, or PySCF. Timings are observed
single runs; campaign timing excludes interpreter startup and does not include
the separate controls, molecular audit, or development time.

The excitation argument also proves a size-independent error ceiling. Its
default value is a loose 2.29046 model units; the direct envelope certificate
gives the much tighter tabled width. Gap controls worsened the direct interval
from 0.01509 at gap 4 to 0.07728, 0.16822, and 0.22205 at gaps 1, 1/4, and
1/20. Doubling hybridization increased it to 0.05920; setting hybridization to
zero reduced it below 10⁻¹⁰. The construction is not uniformly accurate as
the gap closes or the retained coupling grows.

## Independent checks and countercontrols

Eighteen exact tests check the collective CAR normalization, fixed-number
spectator minima, invalid certificate refusals, and resource limits. A direct
rational embedding of all 70 four-particle states of the eight-mode retained
system into twelve physical modes checks its full form matrices, including
bright/dark signs, charging terms, metric factors, and both envelopes.

Six independently assembled full-Hamiltonian controls use at most twelve
physical modes and a sector dimension of 924. Exact operator differences are
diagonal and nonnegative. Numerical ground energies of the full envelopes
agree with the compressed envelopes, and every original ground energy falls
inside its accepted interval. These numerical checks are validation only.
The full-sector routines took 1.679 seconds and 1,523,733 word/state actions;
their additional certificate proposal/replay cost is separate. They refuse
larger inputs before constructing the full Hamiltonian.

A dispersed bath with rank-one coupling has full eight-dimensional Krylov
rank in an exact countercontrol. A small coupling matrix by itself therefore
does not justify exact dark elimination. The operator sandwich is essential
when the original bath is dispersed.

## What happens on the molecular Hamiltonians

The audit examines all 70 H4 and 495 H6 choices of four existing spin orbitals.
It retains all quadratic terms, all active-only quartic terms, and fitted
collective bath and active–bath density interactions. Coefficient medians
minimize the two uniform-density fitting costs. No orbital optimization is
performed, and arbitrary retained quadratic terms are not asserted to fit
the two-band model.

| Frozen molecular fixture | Selected active spin orbitals | Exact lower bound on selected residual norm | Sufficient upper bound on residual norm |
|---|---|---:|---:|
| H4, 8 modes / 4 electrons | 0, 2, 4, 6 | 0.22803364 Ha | 5.51609104 Ha |
| H6, 12 modes / 6 electrons | 4, 6, 8, 10 | 0.28745872 Ha | 20.44544409 Ha |

The lower bounds come from exact off-diagonal matrix elements in the physical
fixed-number sectors. They remain valid after adding a scalar to the
residual. The upper bounds use Hermitian monomial pairs and are conservative.
They would consume twice those upper bounds in a norm-based interval
transfer, far beyond a 1.6 mHa budget. **These are residual norm bounds, not
molecular energy intervals or measured ground-energy errors.**

The selected remainders are therefore not tiny operator perturbations, and
none of the scanned subsets gives a useful sufficient norm-transfer budget.
This does not rule out other orbital choices after rotation, spatially
extended interaction patterns, fixed-number identities, or collective
positivity rules that remain useful even when the residual norm is large.

The quartic tensors also have exact one-leg ranks 8 and 12. This only rules
out representing the same full-Fock tensor on four orbitals by orbital
rotation alone. The twelve-mode compact control itself has full rank before
its fixed-number charging identity is used. Rank is thus not a rejection
test for the broader compression hypothesis. The complete molecular audit,
including exact matrix-element witness discovery, took 10.305 seconds and
explicitly visited the small molecular sectors for those witnesses.

## Decision for the next pass

Keep the original objective. This pass establishes the mechanism when the
collective structure is supplied, and identifies the missing molecular
obligation: collectively handling the nonuniform two-electron interactions.
Further bath enlargement or another local-moment ranking benchmark would not
resolve that obligation.

The next bounded experiment should extract spatially extended **two-electron
operator patterns from H6 itself**, retain their important couplings, and seek
a short operator inequality for the remaining interaction. Candidate tensor
factorizations must reconstruct the fermionic Hamiltonian through exact CAR
algebra, including contraction terms. Success requires a useful certified
remainder together with small pattern description, construction, solve, and
replay costs. A good low-rank fit or a small accepting matrix alone does not
meet that gate. If the proposed grammar fails, record an exact obstruction
for that grammar rather than declaring general many-body compression
impossible.

## Antecedents and artifacts

This is not the first collective elimination in the repository: implicit
Dicke–Schur certificates already appear in
[marginal_implicit_certificate_results.md](/Users/aidenlippert/Documents/Spectra/research/marginal_implicit_certificate_results.md).
The new control tests a particle/hole operator envelope and collective charge
elimination with dispersion. The broader impurity-compression precedent is
[Bravyi and Gosset, *Complexity of quantum impurity problems*](https://arxiv.org/abs/1609.00735).
Neither that precedent nor this construction establishes the required
molecular representation. No literature novelty claim is made; the source
audit is [SOURCES.md](/Users/aidenlippert/Documents/Spectra/research/collective_interference_20260913/SOURCES.md).

Executable sources are in this directory. Data and receipts are in
[results/collective_interference_20260913](/Users/aidenlippert/Documents/Spectra/results/collective_interference_20260913).
Reproduce the accepting and validation checks from the repository root:

```sh
python -S -m unittest research.collective_interference_20260913.test_collective -v
python -S -m research.collective_interference_20260913.replay \
  --campaign results/collective_interference_20260913/campaign \
  --out /tmp/collective-fresh-replay.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  -m research.collective_interference_20260913.validation \
  --out /tmp/collective-full-validation.json
python -S -m research.collective_interference_20260913.scope \
  --out /tmp/collective-molecular-scope.json
```

The campaign command is `python -m research.collective_interference_20260913.campaign
--out <new-directory>`; it refuses to overwrite an existing campaign directory.
