# v10 exponential readout

`experiments/v10_exp_readout.py` evaluates a finite sum

\[
 \sum_j (c_{j,r}+i c_{j,i})T^{d_j}\exp((a_j+i b_j)T)
\]

at nonnegative rational `T`. Every input and every intermediate is a
`Fraction`; positive real exponents are rejected. The returned radius uses the
l1 complex norm `|x|+|y|`, so the exact result lies within that radius of the
returned midpoint.

For an exponent `z`, the implementation chooses a power of two `s` with
`|z/s|_1 <= 1/2`, evaluates `exp(z/s)` by exact complex Taylor arithmetic, and
then repeatedly squares its rational disk. If `u` is the final Taylor term,
the omitted tail is bounded by

`(3/2)|u| / (1 - |z/s|_1/(n+1))`.

The factor 3/2 bounds `exp(|z/s|_1)` for `|z/s|_1 <= 1/2`; the denominator is
the geometric bound on later term ratios. Disk multiplication uses
`r = |m1|_1 r2 + |m2|_1 r1 + r1 r2`, which is valid by submultiplicativity.
Polynomial and sum operations add their propagated radii. Equal exponents are
cached. The requested tolerance is split over unique exponents and weighted
by each polynomial coefficient norm, so the reported total radius is bounded
by tolerance (subject to the exact rational arithmetic and limits).

Fail-closed limits are 50,000 terms, degree 32, 8,192-bit rational limbs, and
100,000 Taylor/scaling iterations. Costs report input terms, unique exponents,
Taylor terms, squarings, and cache hits.
