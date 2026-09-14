# Commutant-obstruction audit

Independent tests pass for the fraction-free inertia routine, integer Gram
reconstruction, and obstruction outcomes. Bareiss/LDL signs match dense
eigenvalue signs on indefinite integer controls; a zero leading pivot is
correctly rejected rather than guessed. For the hidden M4 fixture,
`integer_gram/denominator**2` agrees with a direct rational commutator Gram,
and the threshold obstruction is exactly zero. For the molecular H4 fixture,
the exact replay gives inertia `(positive, negative)=(33,3)`, squared lower
bound `1/8000`, and norm lower `0.011180339887`, above `0.0016`.

Additional exact controls now cover unit-lower-triangular congruences with
known diagonal inertia, the antisymmetric CAR generators, and the
number-ideal core. The core is idempotent, removes a deliberately added
quartic number-ideal lift, preserves diagonal-density structure, and commutes
with a rational mixed orbital rotation on the hidden fixture. The optional
complex-unitary H4 result remains total negative inertia 4, squared lower
`1/10000`, and norm lower `0.01`.

The molecular run initially exposed a portability issue from decimal hashing
of very large Bareiss pivots. The implementation now hashes hexadecimal
integers, and the independent test no longer needs to alter Python's global
integer-string limit.

The theorem's norm is explicitly the Frobenius norm of the quartic
two-particle coefficient representation. It is not a many-body operator norm,
ground-energy error, or coefficient-ℓ1 norm unless a separate conversion
inequality is supplied. The result applies to real orthogonal rotations and
the density-density quartic target family. The antisymmetric CAR generators and
normalized metric were independently checked as well. With the optional
complex-unitary extension, molecular H4 has total negative inertia 4, squared
lower `1/10000`, and norm lower `0.01`. The result still does not address
higher-body terms, many-body operator norm, or broader SOS certificates.
