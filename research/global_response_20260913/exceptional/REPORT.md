# Exceptional compact-sector diagnostic

This pass tested a concrete compact occupation subspace independently of the
previous two-orbital joint bound.  The exceptional sector is defined by
occupation at least three in the last three spatial orbitals.  Its lower bound
is constructed from density-square tangents and a small orbital PSD witness;
no determinant basis or many-body matrix is built.

The result is a precise negative diagnostic.  On frozen H6, the replayed lower
endpoint is approximately **−6.607104 Ha**, while the target is approximately
**−6.333059 Ha**.  On H8, it is approximately **−9.544515 Ha**, also far below
the target.  The occupation pattern therefore cannot certify the exceptional
sector needed for a global interval.  This is a failure of this structured
bound, not a claim that the physical sector has a negative spectral gap.

The test does not reuse the earlier 0.040238-Ha joint gap: it changes both the
number of last orbitals and the required occupation.  It also does not claim a
complement or coupling certificate.  Since the exceptional-sector lower bound
already fails, pursuing those additional pieces would not produce a valid full
bound under this pattern.

The exact checker is `replay.py`; focused tests verify both frozen H6/H8 replay,
zero many-body enumeration, and rejection of a tampered sector partition.
