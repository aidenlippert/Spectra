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
then repeatedly squares its rational disk. If `u` is Taylor term n, put
`r=|z/s|_1/(n+1)`. Submultiplicativity bounds the next term by `|u|_1 r`,
and all subsequent ratios by r. Thus the omitted tail is bounded by

`|u|_1 r/(1-r)`.

No bound `exp(1/2)<=3/2` is used; that claim in the original agent draft was
false and has been removed. Disk multiplication uses
`r = |m1|_1 r2 + |m2|_1 r1 + r1 r2`, which is valid by submultiplicativity.
Each rational midpoint is rounded to a dyadic grid with its exact l1 rounding
distance added to the radius. The radius is rounded upwards. After squaring,
the final disk is explicitly tested against the requested exponent tolerance;
bounded refinement tightens the intermediate tolerance/grid if needed.

Polynomial coefficients are exactly aggregated at T for each (word,exponent)
before exponential evaluation. Let w_e be the sum of these complex l1 norms
over words, and N the number of nonzero exponent groups. Compute each shared
exponential with error <= tolerance/(N max(1,w_e)). Then the summed word
error is <= sum_e w_e tolerance/(N max(1,w_e)) <= tolerance. No assumptions
about cancellation of readout errors or exponent independence enter this sum.

Every eigenoperator word has operator norm1. Summed coefficient l1 error
therefore bounds operator-norm error of the full sum. Since the exact checked
expansion is Hermitian, replacing the approximate operator A by (A+A†)/2
does not enlarge that error. Exact eigenbasis-to-Pauli conversion followed by
taking real Pauli coefficients implements this symmetrization. Truncation and
readout bounds are added at the same operator interface.

Fail-closed limits are 512 words, 50000 total terms (checked while consuming
iterators), degree32, rational bit lengths8192, 128 Taylor terms per attempt,
16 scalings, 16 precision refinements and grid precision4096. Exact rational
inputs are required; floats, strings and bool degrees are rejected. Costs
report input terms, unique exponents, coefficient aggregation, Taylor terms,
squarings, grid rounding and precision refinement. The incorrect agent module
and weak initial tests are retained in `invalidated/`; root rewrote the module
and verified independent rational-series enclosures and refusal behavior.
