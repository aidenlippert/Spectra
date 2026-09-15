# Bounded physical-model branch: the CH2 singlet–triplet separation

The defined observable is E(1A1)−E(3B1). The preserved minimal-basis,
frozen-core calculation gives an exact interval at one fixed geometry.
Its narrow solver interval does not certify an adiabatic electronic or measured
vibrational separation. The higher-basis controls use separately optimized
CASSCF(6,6) singlet and triplet geometries; SC-NEVPT2 and CCSD(T) are single
points on those geometries. CCSD(T) correlates all eight electrons, whereas the
CASSCF model treats six active electrons with one inactive doubly occupied orbital.

The [primary photoelectron study](https://doi.org/10.1063/1.449746) reports an
observed splitting of 9.00±0.09 kcal/mol and infers an electronic separation
of 8.7±0.5 kcal/mol after considering vibrational effects. These are different
comparison quantities. This campaign does not calculate zero-point corrections.

The machine-readable ladder reconstructs the exact preserved interval using
[L_S−U_T, U_S−L_T], records the numerical controls and binds every input by its
hash. Its kcal/mol values are display conversions; the rational hartree endpoints
remain the exact mathematical result. The exact model solver width is roughly
0.013 kcal/mol, while the model discrepancy is much larger.

Comparing the double- and triple-zeta rows changes both the basis and its
optimized geometry. That change is a sensitivity, not a rigorous error estimate.
Differences between CASSCF, NEVPT2 and CCSD(T) similarly are not confidence bounds.
No combination of those differences is promoted to a certified physical interval.

A useful subsequent physical calculation would establish compatible geometry
and nuclear corrections for the selected electronic method, then test convergence
of the model against the electronic observable. It is secondary to the interacting
certificate campaign. No experimental recommendation or advantage is claimed here.
