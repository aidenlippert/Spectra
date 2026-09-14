# Hidden density-basis audit

Independent controls in `test_hidden_density_basis.py` pass. For every
generator on a four-mode hidden fixture, the implementation's CAR commutator
column agrees with a direct symbolic commutator using the derivation rule.
The rotated Hamiltonian agrees with explicit word-by-word substitution under
two overlapping rational Givens rotations, including quadratic and quartic
terms; inverse transport recovers the original Hamiltonian. M4 H-only discovery
produces an exactly orthogonal transport certificate whose inner ratio-chain
replay has zero width. Bad rotations, sector mismatch, Hamiltonian mismatch,
and corrupted inner upper metadata are rejected.

The commutator shortcut is therefore validated for the tested normal-ordered
quartic CAR words, including the ascending annihilator convention. The exact
replay checks rational orthogonality, transformed-H identity, and the inner
certificate; numerical nullspaces and eigenspaces remain discovery proposals,
as recorded by the implementation. For the M4 hidden fixture, an exact
modular minor gives rank at least 6 and hence nullity at most 4, matching the
observed numerical nullity 4. This certifies the dimension bound, not the
floating nullspace basis itself.

Scope is narrow and correctly stated: this transports a known solvable ratio
family through a real rational orthogonal orbital basis. It does not discover
generic hidden structure, prove a general commutant theorem, or address
molecular Coulomb Hamiltonians. The M4/M6/M8 size runs demonstrate the cost of
this recognition-and-transport path conditional on the family.
