# Exact occupation degree of the valence projector

On m=2p spatial sites with Nα=Nβ=p, let D count doubly occupied sites. The indicator P of the singly occupied valence space has minimal occupation-polynomial degree m, even when equality is required only on this fixed-number Boolean slice.

An upper representation is

\[
P=\prod_{j=1}^{p}(1-D/j).
\]

Every physical occupation has D in {0,…,p}, and every such value occurs. Since D has occupation degree two, this gives degree 2p=m.

For the lower bound, average any representing polynomial over permutations of the spatial sites. This preserves its degree and the represented function. A monomial with a alpha occupations and b beta occupations, including c common sites, averages to a polynomial in D of degree at most min(a,b). Explicitly, with (x)_k the falling factorial, its average at D is

\[
\frac{(D)_c}{(m)_{a+b-c}}
\sum_{j=0}^{\min(a-c,b-c)}(-1)^j
{a-c\choose j}{b-c\choose j}j!
(D-c)_j(p-c-j)_{a-c-j}(p-c-j)_{b-c-j}.
\]

The factors containing p are constants. The D degree is at most c+min(a−c,b−c)=min(a,b), which is at most floor((a+b)/2). Monomials with a>p or b>p vanish on the slice and can be removed.

Thus an occupation polynomial of degree r representing P would yield a univariate polynomial of degree at most floor(r/2), with value 1 at D=0 and zeros at all p distinct points 1,…,p. It must have degree at least p. Therefore r≥2p=m.

The exact averaging formula and upper representation have also been checked by explicit rational counting at 4, 6 and 8 sites: 42, 120 and 275 averaging identities respectively, over 36, 400 and 4,900 occupations. These checks support the formula; the general result follows from the argument above.

This theorem concerns expanded occupation polynomials for the exact projector. It does not exclude a factorized product, a compact tensor operator, a vanishing metric, or other spectral certificates. In particular, the metric condition v(0)=0 removes the target projector from the weighted-row numerator altogether. [Construction and measured limits](marginal_joint_metric.md).
