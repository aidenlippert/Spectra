# Stoquastic chain certificate audit

Independent dense controls in `test_stoquastic_chain_certificate.py` pass for
periodic chains of 4 and 5 sites. The local-energy dynamic program matches
brute-force configuration extrema, and the exact `[lower, upper]` interval
contains the dense Hamiltonian ground eigenvalue for the parent family at
three rational `q` values. A perturbed diagonal table and nonuniform positive
flip rates are also bracketed. Independent local-term matrices have a
nonzero commutator, so the control is not merely a commuting diagonal model.
Negative flip rates are rejected.

The certificate proves a restricted stoquastic positive-amplitude/Jastrow
family. It does not prove a generic fermionic or spin-chain ground-state
algorithm. The optimizer in `discover()` is a floating `minimize_scalar`
heuristic followed by rationalization; its success flag and local objective
do not prove the best `q` or a zero-width result for the family. The verifier
itself checks the exact conditional 2x2 PSD identities and DP extrema.

The local family has constant-size tables and a cyclic DP, so its replay work
is linear in chain length for fixed local factor width and rational bit length.
That is a family-specific structural result. It should not be generalized to
unrestricted interacting fermions, arbitrary sign problems, or a polynomial
ground-state theorem. The large zero-width ladder runs demonstrate this
restricted construction's behavior only.
