# Fixed-number CAR kernel and the number-ideal lift

Let `P` be a number-conserving, normally ordered CAR polynomial of operator
degree at most six (body rank at most three), and restrict it to the `N`-
particle sector of `M` modes. The useful symbolic question is whether a zero
sector operator can be represented as

`P = (Nhat-N) X + R`,

with `X` number-conserving of degree at most four and `R=0` after CAR reduction.

For the generic range `3 <= N <= M-3`, this is the natural degree-bounded
number ideal statement: expand `P` in the normal-ordered body basis
`a†_{i1}...a†_{ik}a_{j_k}...a_{j1}`, `k=0,1,2,3`; multiplication by
`Nhat-N` maps body ranks `0,1,2` in `X` to ranks at most `1,2,3` after
contractions. Solve the resulting exact linear coefficient map. Its matrix is
an incidence/contraction map on subsets, with dimensions governed by
`binom(M,k)^2` (and symmetry blocks), rather than by the sector dimension
`binom(M,N)` or the full Fock dimension `2^M`.

There is an important qualification. This is a statement about the kernel of
the *representation map on the N-sector*, after quotienting the CAR relations;
it is not an assertion that every arbitrary six-degree word has a unique lift.
A lift can be nonunique if multiplication by `Nhat-N` has a kernel within the admitted multiplier space; this is not assumed for M=6,N=3. The implemented numerical map has full reported column rank in that case, without using this observation as a proof. A practical verifier constructs the rational matrix of
the coefficient map, computes rank and residual by exact Gaussian elimination,
and reports a basis of the cokernel. No floating least-squares result is a
certificate.

The safe rank criterion is controlled by both particles and holes. A k-body
operator basis can fail to be faithfully represented when `N<k` or
`M-N<k`; in those regimes top-rank tensors can vanish on the sector or collapse
through finite-size identities. For `N=3`, rank-three terms sit at the
threshold, so the generic claim must be checked explicitly for each M and
symmetry block. At `N=0,1,2` (or fewer than three holes), one must lower the
body basis and include the resulting boundary identities.

For the compiler, this gives an exact symbolic lift routine: canonicalize words,
collect body-rank coefficients, build the contraction/incidence matrix for
`(Nhat-N)X`, solve for a rational multiplier, and retain any nonzero cokernel component in an explicitly norm-bounded residual. Exact equality requires zero residual; a safe approximate lower bound does not.
The routine can certify a proposed SOS identity's number-sector equality while
keeping the positivity part separate. It does not by itself prove a low-rank
certificate exists or make the marginal cone easy.
