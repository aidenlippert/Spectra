# Obstructions to a compact Schur-response theorem

This note is a stress test for the proposed response

\[
 \Sigma(E)=B^*(C-E)^{-1}B,
 \qquad C=QHQ,\;B=QHP,
\]

where the target is to construct certified upper and lower response bounds
without first solving the full eliminated sector.  The examples below are
finite-dimensional and therefore exactly checkable.  They do not prove that
Spectra's construction cannot work on its molecular class; they identify
hypotheses that are insufficient on their own.

## A gapped eliminated sector can have an extensive response

Let `dim(P)=1`, let `Q` have orthonormal basis `q_1,...,q_N`, and define

\[
 C=\operatorname{diag}(\Delta_1,\ldots,\Delta_N),\qquad
 Bp=\sum_{j=1}^N g_jq_j,
\]

with `\Delta_j\ge \Delta>0` and `E\le 0`.  The eliminated sector has a
uniform gap at least `\Delta` above the energy window, but

\[
 \Sigma(E)=\sum_{j=1}^N\frac{|g_j|^2}{\Delta_j-E}.
\]

If the `\Delta_j` are distinct and every `g_j` is nonzero, this is a rational
function with `N` distinct poles.  Its exact specification requires all
weights and pole locations in general.  This is an exact-representation
obstruction, not a fixed-tolerance approximation lower bound: away from the
poles, polynomial or rational approximation can sometimes have degree
independent of `N`.  A gap controls the norm,

\[
 0\preceq\Sigma(E)\preceq \frac{\|B\|^2}{\Delta-E}I,
\]

but does not control the number of poles or the number of moments needed for a
uniformly certified approximation.  In particular, repeated elimination can
replace one small retained block by a memory kernel with an extensive pole
set, even though every eliminated excitation is individually expensive.

The same obstruction can be made matrix-valued.  Take `dim(P)=N`, `C` as
above, and choose `B` with dense columns (for example a scaled orthogonal
matrix times a diagonal of nonzero couplings).  Then

\[
 \Sigma(E)=B^*\operatorname{diag}((\Delta_j-E)^{-1})B
\]

is generically full rank and dense for every `E`.  A spectral gap alone does
not imply a small number of collective boundary channels, a sparse induced
operator, or a low-rank response.

## Locality and area-law language do not repair the implication

The star construction above is a local interaction on a graph whose retained
site is adjacent to each eliminated site; graph locality by itself permits a
boundary with `N` channels.  A small geometric boundary is therefore not the
same assumption as a small *operator channel rank*.  More generally, higher-
dimensional area laws do not by themselves imply an efficient classical
description, and PEPS contraction is hard in the worst case.  These are
complexity cautions, not claims about every physical instance:

* Schuch, Verstraete, Cirac, and Wolf, *Computational complexity of projected
  entangled pair states*, arXiv:quant-ph/0611050.
* Ge and Eisert, *Area laws and efficient descriptions of quantum states*,
  arXiv:1411.2995.

Likewise, a low bond dimension or a small Schmidt tail for one prepared state
does not bound the response of an arbitrary eliminated-sector excitation.  A
state-dependent response theorem must state the observable, the initial state,
and the weighted omitted-state error explicitly.

## Exact retained-weight theorem

Let `H\psi=E_0\psi`, let `P+Q=I` be orthogonal projectors, and write
`u=P\psi`, `v=Q\psi`, `p=\|u\|^2`, and `C=QHQ`, with `0<p<1`.  The projected
eigenvalue equations imply the exact identity

\[
 \langle v,(C-E_0)v\rangle
 =\langle u,(PHP-E_0)u\rangle.
\]

Therefore,

\[
 \lambda_{\min}(C)-E_0
 \leq \frac{\langle u,(PHP-E_0)u\rangle}{1-p}.
\]

For a half-filled Hubbard model with `P` the global no-doublon projector,
`PHP=0`, hence `\lambda_{\min}(C)-E_0\le -E_0p/(1-p)` when `E_0<0`.

For `L` independent two-site Hubbard dimers with `U=8,t=1`, the exact dimer
ground energy and no-doublon weight are

\[
 E_d=4-2\sqrt5,\qquad w=\frac{1+2/\sqrt5}{2}=\frac{5+2\sqrt5}{10}.
\]

The product ground state has `E_0=LE_d` and `p=w^L`; thus a fixed local
charge gap does not prevent the global bare no-doublon weight from decaying
exponentially.  Moreover, from
`(C-E_0)^{-1}QHP\,u=-v`,

\[
 \|(C-E_0)^{-1}QHP\|\ge\sqrt{(1-w^L)/w^L}.
\]

This only obstructs one global bare-sector elimination.  It does not rule out
dressed or local adaptive projectors; those must control response channels
compositionally.

## What this theorem does not prove

The exponentially small `p` is not by itself an exponential lower bound on
the cost of approximating the response at a fixed absolute energy tolerance.
If the evaluation point is separated from the eliminated spectrum by
`\varepsilon>0`, then `\|(C-E)^{-1}\|\le 1/\varepsilon`, and polynomial or
rational approximation may have degree independent of the eliminated-sector
dimension for suitable spectral geometries.  Also, with
`X(E)=(C-E)^{-1}B`,

\[
 \frac{d}{dE}\bigl[B^*(C-E)^{-1}B\bigr]=X(E)^*X(E)\succeq0.
\]

Thus the Feshbach residual slope contains `I+X^*X`; a large response slope
can convert response error into a smaller eigenvalue error.  The dimer result
identifies a failure mode for a bare global overlap argument, not a complexity
lower bound.  A genuine obstruction would need to prove large weighted rank,
verification work, or composition growth for every certified construction in
the declared model class.

## What would escape the counterexample

The needed assumption is stronger than a gap.  A useful, independently
checkable condition could be a **weighted boundary-response complexity bound**:

\[
 \mu_B=\sum_j |g_j|^2\,\delta_{\Delta_j},\qquad
 \sup_{E\in I}\left\|
 B^*(C-E)^{-1}B-\widetilde\Sigma(E)
 \right\|\le \varepsilon,
\]

where `\widetilde\Sigma` has at most `r(\varepsilon,I)` certified poles or
channels, and both `r` and the certificate construction are bounded from
input data before diagonalizing `C`.  Equivalent usable conditions include a
certified decay bound on the Stieltjes moments

\[
 m_k=B^*C^{-(k+1)}B,
\]

or a directly computed bound on the numerical rank of the weighted Krylov
space

\[
 \mathcal K_r(C,B)=\operatorname{span}\{B,CB,\ldots,C^{r-1}B\},
\]

with the residual of the rational/Krylov approximant bounded uniformly on
`I`.  The rank and residual must be obtained using local operator algebra or a
certified iterative procedure; fitting after constructing the full response
does not establish cheap discovery.

For a composable theorem, the condition must survive elimination: each new
boundary response needs a certified channel/rank bound and an additive error
rule, while the induced retained Hamiltonian remains representable in the
same class.  A practical campaign should therefore measure, before claiming
compression: (i) weighted channel rank versus tolerance, (ii) uniform
resolvent residual over the energy interval, (iii) boundary size and support,
(iv) growth after a second elimination, and (v) cold construction cost and
peak memory.

## Complexity claim to avoid

“The eliminated sector is gapped,” “the state obeys an area law,” or “the
response has a small file” is not a complexity theorem.  A valid theorem must
give a physically recognizable input class and prove that the response rank,
approximation error, verification work, and repeated-composition cost are
bounded for that class.  Worst-case fixed-basis electronic-structure
Hamiltonians are QMA-complete (O'Gorman et al., arXiv:2103.08215), so a
universal cheap construction would require a substantially stronger restriction
than generic locality or a gap.

## Constructive conclusion

The counterexample narrows the breakthrough target. The result to seek is not
“Schur elimination is compact for gapped systems.” It is a theorem identifying
a physically important class with a pre-solution, certifiable bound on the
weighted boundary spectral complexity, together with a response constructor
and a composition theorem. If the measured H8/H10 or coupled-molecule cases
show small weighted rank and stable rank under a second elimination, they are
evidence for that hypothesis; they are not yet a proof of it.
