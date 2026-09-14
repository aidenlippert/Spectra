# CH₂: exact spin-gap certificate for a specified small model

The pilot now has independently replayed **S=0 and S=1 energy intervals**.
For the declared fixed-geometry, frozen-core STO-3G Hamiltonian,

    E_T − E_S ∈ [−25.49213611, −25.47958594] kcal/mol (approximate conversion).

The exact interval is stored in rational hartrees in
`results/ch2_validation_20260913/exact_replay.json`; its width is
**0.0199999583 mHa**. The triplet is lower in this model. This demonstrates
solver and state-identification validation, not accurate prediction of the
experimental splitting or an advantage over existing chemistry solvers.

## Model and states

The fixed coordinates, in Å, are C=(0,0,0), H=(0,0,1.117),
H=(0,1.047,−0.389). Neutral CH₂ has eight electrons. STO-3G produces seven
spatial orbitals here. One RHF core orbital is frozen, leaving **six active
electrons in all six remaining spatial orbitals**, CAS(6,6).

The importer uses PySCF 2.14.0 `CASCI.get_h1eff/get_h2eff`. This includes the
frozen-core Coulomb and exchange potential in the active one-body operator,
and the frozen-core energy plus nuclear repulsion in the common constant.
That treatment follows the
[PySCF CASCI source](https://pyscf.org/_modules/pyscf/mcscf/casci.html).
The RHF total is −38.3687284528834 Ha and the common constant
−27.9080109317997 Ha. A focused test reproduces the RHF expectation from the
exported Hamiltonian within the rounding envelope.

Spatial one- and two-electron coefficients are rounded symmetrically at
denominator 10¹² before spin expansion. This preserves the spin symmetry.
The importer computes a conservative coefficient-norm rounding envelope of
about 4.086e−10 Ha against its floating-point integral inputs. It does not
provide interval enclosures for integral evaluation, orbitals, or the physical
model. The certificates below concern the exported rational Hamiltonian.

Both states use this same geometry and orbital definition. No state-specific
geometry optimization, zero-point energy, relativistic, thermal or solvent
correction is included.

## Exact state and energy proof

The accepting path explicitly constructs the complete M_S=0 block
(400 determinants) and M_S=1 block (225 determinants). It constructs S²
from the exact fermionic spin operators and checks [H,S²]=0 in each block.

For S=s, a rational factor proves

    H + 20 Ha [S²−s(s+1)I] − L_s I ≥ 0

on the complete M_S=s block. On the target S=s eigenspace the penalty is
zero, which proves the lower bound on that total-spin sector. The upper
trial is projected by an exact polynomial in S², and replay verifies
S²v=s(s+1)v before accepting its Rayleigh quotient. M_S alone is never
used as a pure-spin certificate. See also
[PySCF FCI documentation](https://pyscf.org/user/ci.html) for the distinction
between spin projection settings and spin-adapted solvers.

Each absolute energy interval is approximately 0.01 mHa wide. The gap
combines them as [L_T−U_S, U_T−L_S]. Exact replay imports no NumPy, SciPy,
PySCF or CVXPY. Five focused tests pass, including incorrect electron/core
accounting, altered factors and spin-contaminated upper-witness refusals.

The independent two-electron/two-orbital state control checks a four-state
M_S=0 space with **three singlets and one triplet component**. It is an
algebra test, not a CH₂ physical model.

## Existing solver and experimental comparisons

PySCF's numerical FCI controls for the same floating-input model give
S=0: −38.429076457852574 Ha and S=1: −38.469690762594965 Ha, with measured
S² near 0 and 2 respectively. Their energy difference agrees with the narrow
certified interval. The two numerical FCI calls take about 0.023 s and
0.006 s in this one run. This is a baseline comparison, not evidence of
Spectra being faster or more accurate on the model.

Leopold, Murray, Miller and Lineberger report a CH₂ singlet–triplet splitting
of **9.00±0.09 kcal/mol**, with the triplet ground state. Their abstract also
reports an isotope-based zero-point contribution of 0.27±0.40 kcal/mol and
an inferred electronic separation T_e of 8.7±0.5 kcal/mol.
[Primary study abstract, University of Minnesota](https://experts.umn.edu/en/publications/methylene-a-study-of-the-xsup3supbsub1sub-and-%C3%A3sup1supasub1sub-st/).

In this report's E_T−E_S convention the observed splitting is negative.
The present roughly −25.49 kcal/mol fixed-geometry minimal-model result
does **not** establish physical accuracy. It also is not the same observable
as a relaxed, zero-point-corrected splitting, so simply subtracting the two
numbers would not give a rigorously assessed model error. Basis, frozen-core,
geometry and nuclear-motion errors require separate work.

## Costs, invalidated attempts and replay

Corrected fixture generation, RHF and numerical FCI controls took 0.143 s
after chemistry imports. Certificate construction plus its initial exact
checks took 2.199 s; a separate exact replay took 1.064 s. These are measured
stage times, with complete determinant matrices explicitly charged. The
full accounting and limitations are in `cost_ledger.json` in the result folder.

Earlier proposals had wrong active-electron/orbital accounting, incorrect
spin labels, or omitted frozen-core one-body contributions. They are
invalidated and archived under `invalidated_core_omission` and
`invalidated_4e_proposal`; none contributes to the accepted result. Their
reported timings were incomplete and are not silently counted as zero.
The old `run_active.py` entry point refuses to produce those outputs.

Read-only exact replay of this pilot together with the main experiment:

```sh
cd /Users/aidenlippert/Documents/Spectra
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.composable_response_20260913.verify
```

This pilot should remain secondary to the retained-sector response work.
Its useful next physical step is a properly matched geometry/basis study,
with the same exact spin identification and a separate model-error assessment.
