# Physical-validation pilot: gas-phase carbene spin gaps

Start with methylene, CH2. The decision is which of its lowest singlet and
triplet states is lower, and the quantitative target is their adiabatic
zero-vibrational-level separation. This is a proposed retrospective pilot;
no Spectra chemical prediction or laboratory advantage is established here.

Use one sign convention throughout: Δ_ST = E_T − E_S. A positive value
means that the singlet is lower. Electronic minima, vibrational origins,
vertical excitations and finite-temperature free energies are distinct
observables and must occupy separate columns.

## Reference records checked on September 13, 2026

| Case | Role | Reference for Δ_ST | What the reference supports |
|---|---|---|---|
| CH2 | First development case | Approximately −9.0 kcal/mol | NIST's photoelectron compilation gives a 9.0 kcal/mol splitting and cites Leopold et al. (1985). It does not supply a gap uncertainty in that entry. Do not mistake the adjacent electron-affinity uncertainty for an uncertainty on the spin gap. [NIST Methylene](https://webbook.nist.gov/cgi/cbook.cgi?ID=C2465567&Mask=20) |
| CF2 | Proposed transfer case after freezing CH2 modeling | +54 ± 3 kcal/mol | Schwartz et al. report resolved photoelectron transitions to the two states, with the singlet lower. This is a historical experimental reference with limited precision. [Schwartz et al., 1999](https://doi.org/10.1021/jp992214c) |
| CCl2 | Optional transfer and reference-quality control | Estimated +0.9 ± 0.2 eV, approximately +20.8 ± 4.6 kcal/mol | Wren et al. corrected contaminated earlier spectra. The triplet origin was not directly observed, so this is an estimated separation rather than a precision band-origin measurement. [Wren et al., 2009, primary abstract](https://pubmed.ncbi.nlm.nih.gov/19492128/) |

**The old CCl2 value of 3 ± 3 kcal/mol is superseded.** The 2009 authors
explicitly identify dihalomethyl-anion contamination and reject their earlier
interpretation. CCl2 must therefore not be advertised as a measured
near-degeneracy or scored against the 1999 value. This is a concrete example
of why a precise solver must also audit its physical reference.
[Wren et al., 2009](https://doi.org/10.1039/b822690c)

The two halocarbenes are candidates for later transfer tests, not a frozen
sub-kcal experimental benchmark suite. Their reference uncertainties cannot
validate a sub-kcal improvement. A prediction matching those records would
initially support broad state ordering and consistency only.

## A concrete first comparison

For CH2, compare Spectra and a fixed canonical CCSD(T) calculation using the
same Hamiltonian, geometry and frozen-core convention. Canonical CCSD(T) is
a proposed baseline, not certified truth. A published study of methylene
and its adducts demonstrates that spin gaps and geometry-dependent shifts
are accessible to established coupled-cluster workflows; Spectra must show
what its guarantee adds. [Altun et al., local energy decomposition study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6728066/)

Before producing a prediction, record the following in a calculation manifest:

1. Both state-specific geometries, basis, active orbitals, frozen core,
   integrals, constants, and their hashes. Use total-spin S=0 and S=1 sectors;
   fixing spin projection M_S alone does not isolate a singlet. The present
   H6 total-particle checker does not yet supply these two spin-sector proofs.
2. A small active-space version on which exact diagonalization can check the
   spin-sector implementation. Its role is solver verification. Separately
   assess basis and excluded-correlation errors before comparing with matter.
3. Separate electronic gaps Δ_e and zero-point corrections Δ_ZPE, with
   Δ_0 = Δ_e + Δ_ZPE. Record the experimental isotopologue and energy-zero
   convention; do not insert a vertical gap into this comparison.
4. Spectra's two exact energy intervals, the matched baseline, all construction
   and replay costs, and the origin of every upper witness. Reference FCI
   discovery remains part of the accounting when used.

For each supplied Hamiltonian, interval subtraction is exact:

    E_S in [L_S,U_S], E_T in [L_T,U_T]
    E_T − E_S in [L_T−U_S, U_T−L_S].

The difference width is the sum of the two absolute widths unless an
additional joint certificate justifies a tighter result. Solver certificates
do not automatically bound model, geometry or zero-point errors.

## Proposed acceptance criteria

Require a solver-only CH2 gap interval no wider than **0.5 kcal/mol**, with
its sign resolved. Use **1.0 kcal/mol** as the initial total physical
comparison tolerance. These are proposed decision tolerances, not measured
error bars. Model and reference-convention uncertainty must be assessed
independently; if their combined contribution cannot fit the remaining
budget, mark the physical comparison inconclusive. First retrieve and audit
the primary experimental gap uncertainty rather than assigning the NIST
entry an invented error bar.

Record whether Spectra delivers a certified answer at competitive cost or
resolves a question that the baseline leaves uncertain. Reproducing the
known CH2 ordering alone is a sanity check. Passing one development molecule
does not demonstrate predictive advantage.

After that first comparison, freeze the procedure before any halocarbene
solver results. Preserve unsuccessful cases. CCl2 remains optional until its
model and reference limitations suit the intended decision. This pilot uses
existing data; prospective laboratory tests and commercial claims remain
later milestones.
