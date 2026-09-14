# Sparse response lift upper bound

This experiment constructs

`psi(alpha) = (I - alpha Q H) |HF>`

where `Q` is double occupation of the highest spatial orbital. The closed-shell
Hartree–Fock determinant occupies the first `particles/2` spatial orbitals in
the interleaved spin ordering, so `Q|HF>=0`. Applying the decoded Hamiltonian
directly to that one bitstring and retaining only Q-supported outputs produces
the finite excitation vector `chi=QH|HF>`. No fixed-particle determinant basis
or inherited CI amplitudes are used.

The exact Rayleigh quotient is

`(e0 - 2 alpha ||chi||^2 + alpha^2 <chi,H chi>) /
 (1 + alpha^2 ||chi||^2)`.

Alpha is selected by an auditable rational grid search with denominator 1000
over [-10,10]. All moments and the final upper are exact rational values.

Results (including fresh H6 geometries at 1.6 Å and 1.73 Å):

| fixture | alpha | chi terms | lifted upper |
|---|---:|---:|---:|
| H6 | 29/40 | 5 | -6.139059287947763 Ha |
| H8 | 661/1000 | 8 | -8.991999023351374 Ha |
| fresh H6 1.6 Å | rational grid optimum | 5 | see certificate |
| fresh H6 1.73 Å | rational grid optimum | 5 | see certificate |

The construction takes under a second for all four cases. Exact replay binds
the fixture hash, alpha, sparse support, both moments, numerator, denominator,
and upper endpoint; mutation tests reject altered alpha. Materialized labels
for `H|HF>`, `QH|HF>`, and `H(QH|HF>)` are reported separately, so sparse
excitation support is not misreported as total work. Independent sparse CAR
tests pass 3/3. The result is an upper-bound response component and does not by
itself provide a tight interval or a physical accuracy claim. Its sparse
support is evidence that the response lift avoids importing the previous
enumerated upper amplitudes; the Hamiltonian word action itself remains charged
in the construction.
