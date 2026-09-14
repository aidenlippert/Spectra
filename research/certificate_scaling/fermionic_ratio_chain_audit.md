# Fermionic ratio-chain audit

Independent controls in `test_fermionic_ratio_chain.py` pass for ratios
`r=1,4,1/4` and every particle sector `N=0,...,4` on four modes. They apply
the complete Hamiltonian polynomial directly to the weighted fixed-sector
vector and obtain an exact zero vector, rather than checking only an energy
expectation. Compilation and structured replay return zero width. Adjacent
local terms have a nonzero algebraic commutator, coefficient mutation is
rejected, wrong upper-state ratio is rejected, and malformed mode/sector
parameters are rejected.

The compiler is correctly H-only for its declared ratio family: it recognizes
the full canonical Hamiltonian, derives the rational square ratio, and emits
cubic CAR factors whose production replay gives lower zero. The structured
upper replay uses an exact elementary-symmetric dynamic program for the norm,
with no many-body state enumeration. Its zero-width scaling results are
therefore credible for this solvable exclusion/ferromagnet-like family.

Scope remains narrow. This is a known structured parent family with a
one-parameter ground state and an exact local factorization. It does not imply
a generic fermionic ground-state compiler, molecular applicability, or a new
general representability theorem. The numeric size sweep demonstrates output
and replay scaling conditional on the family; it does not establish discovery
of such structure for arbitrary Hamiltonians.
