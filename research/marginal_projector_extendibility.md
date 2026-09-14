# Quantitative extension constraints from overlapping projectors

Follow-up: [exact hopping filters](marginal_hopping_filter.md) now correlate all
neighboring eight-site blocks and tighten the million-site interval to
[-0.611636,-0.5581139722], with the displayed upper rounded outward. The
projector bounds below remain valid; the uncorrelated upper is an earlier stage.

The four-site compatibility obstruction now gives exact linear inequalities,
including a stronger Hubbard energy certificate. The construction breaks the
previous four-site boundary-correction ceiling using local matrices of dimension
at most20 and projector Gram blocks of dimension at most36.

For the half-filled uniform Hubbard chain at U=4,t=1, the refined certificate gives

\[
E_0^{\rm periodic}/N\ge
-\frac{37473363594937}{60000000000000}
\approx-0.624556059915617.
\]

For an even open chain with N≥8, subtract2/N. At one million sites the lower
energy density is approximately-0.624558059915617. The separately certified
six-site window still gives the stronger lower bound-0.611636 there. This new
result addresses the missing extension constraint with smaller matrices; it
does not replace that best energy bound.

## The geometric inequality

Let φ be the integer-amplitude four-site pure witness in
`results/marginal_graded_hubbard8/weighted_window_family/certificate.json`,
normalized by its squared norm n. Write P_j for its rank-one projector embedded
on sites j through j+3. Although φ has identical left and right three-site
reductions, consecutive copies cannot simultaneously have fidelity1.

For m consecutive windows on m+3 sites, write P_j=V_j V_j†. Each V_j is an
isometry with4^(m−1) columns, obtained by adjoining all exterior occupation
states. If W=[V_1 ... V_m], then

\[
\sum_{j=1}^m P_j=WW^\dagger,
\qquad G_m=W^\dagger W.
\]

The nonzero spectra coincide. Thus an exact Gram inequality
G_m≤B_m I proves ΣP_j≤B_m I on the entire local Fock space. The verifier builds
every sparse isometry column. The supplied witness has fixed spin numbers(2,2),
so different total spin-number sectors are exactly orthogonal and may be checked
separately. For a submitted vector without that symmetry, only total particle
number is used to split the Gram.

| Windows m | Total Gram dimension | Largest checked Gram block | Exact accepted B_m | Fidelity cap B_m/m |
|---:|---:|---:|---:|---:|
| 2 | 8 | 2 | 1.380638468 | 0.690319234 |
| 3 | 48 | 12 | 1.901418483 | 0.633806161 |
| 4 | 256 | 36 | 2.223786409 | 0.55594660225 |

Every number in the last two columns is a rational certificate endpoint, not
an accepted floating-point eigenvalue. The routine clears denominators and
checks B_m n I minus the integer Gram by fraction-free PSD elimination.
It admits only m=2,3,4 and retains the existing matrix budget.

The two-window result has an additional direct check: V_1†V_2 is diagonal, with
entries

\[
\left[
\frac{106854722802607664801035}{895219092079179430700602},
\frac{170377411618491025274633}{447609546039589715350301},
\frac{170377411618491025274633}{447609546039589715350301},
\frac{106854722802607664801035}{895219092079179430700602}
\right].
\]

Its norm is about0.380638467445513, yielding the same two-window bound.
The tests independently embed the two isometries in five-site occupation space
and check the unsplit8×8 Gram inequality exactly.

## From finite overlap to the chain

For a periodic chain longer than m+3 sites, cyclically translate the local
inequality. Each P_j occurs exactly m times, hence

\[
\sum_{j=1}^N P_j\le N\theta_m I,
\qquad \theta_m=B_m/m.
\]

No translation symmetry of the state is assumed. For the average four-site
marginal, this proves Tr(P_φ ρ_average)≤θ_m. At m=4 the excluded fidelity margin
is0.44405339775. In the conventional trace distance D(ρ,σ)=||ρ−σ||_1/2, testing
with the effect P_φ also gives D(ρ_average,P_φ)≥0.44405339775. This excludes a
neighborhood of the old witness, not only the witness itself.

All contiguous projectors are even fermionic operators. Wrapped copies are
defined by conjugating with the actual fermionic mode-translation unitary;
one must retain its signs in an occupation-basis implementation. The proof
uses unitary conjugation and therefore does not assume that a naive tensor
permutation gives the wrapped matrices.

## The energy certificate

Use the four-site centered window with onsite profile[a,6−a,6−a,a] and hopping
profile[b,3−2b,b]. Its cyclic sum equals3 times the centered uniform Hubbard
operator at U=4,t=1. For κ≥0, check

\[
K_4+\kappa P_\phi\ge\ell I
\]

in every local particle and reflection sector. Combining with the projector
inequality gives

\[
\frac{E_0^{\rm periodic}}N\ge\frac{\ell-\kappa\theta_m}{3}.
\]

Explicitly, three times the claimed shifted bulk operator is the sum of
translated positive terms K_4+κP−ℓI, plus κ/m times the translated positive
terms B_m I−Σ_{r=0}^{m−1}P_{j+r}. All projector coefficients cancel. At half
filling the centered bulk shift vanishes. Removing the periodic wrap hopping
costs at most its operator norm2, giving the open-chain correction-2/N.

The refined rational parameters are

\[
a=0.5613,\quad b=0.8171,\quad\kappa=0.4386,\quad
\ell=-1.62983,\quad m=4.
\]

All46 local blocks pass exact positivity and cover all256 local Fock states.
The largest local block is20×20. The energy verifier explicitly requires the
actual submitted projector to preserve spin-number and reflection sectors;
otherwise those separate block checks would be invalid.

The original profile a=0.531373,b=0.75 also passes with κ=0.3157 and
ℓ=-1.724725:

| Constraint | Certified periodic lower density, approximately |
|---|---:|
| Old four-site boundary family: upper ceiling on achievable lower bound | -0.680141558183 |
| Two-projector inequality | -0.647552927391 |
| Three-projector inequality | -0.641605868343 |
| Four-projector inequality | -0.633412447443 |
| Four-projector inequality with refined profile and penalty | -0.624556059916 |
| Separate six-site window lower bound | -0.611634 |

The old ceiling applies even with arbitrary Hermitian even three-site boundary
differences. Therefore the improvement is evidence of a new admissible
extension constraint, beyond further optimization within that old family.

An intermediate numerical proposal used only local half filling and was
discarded: other local particle sectors had lower energy. The final certificate
is accepted solely by the all-Fock exact replay. Numerical search is not an
optimality proof; one exploratory agent also exceeded its requested numerical
evaluation budget, so no bounded-search performance claim is made.

## A stronger extensive physical upper

The previously verified eight-site polynomial witness now has an upper-only
replay and an explicit tiling adapter. Its energy is

\[
u_8=-\frac{120558410244999997006043477431315508062486123707}
{28461740687412294825487447957163629238443780483}.
\]

The adapter validates the original uniform H8 CAR operator, each fixed-particle
boundary vector, and the bounded integer Chebyshev table. It recomputes moments
through order23 and the exact physical Rayleigh quotient. The retained polynomial
witness needs no ground-spin theorem or complement-gap proof to establish its
variational upper bound. It preserves the existing symmetry-oracle limits.

For N=8q+r with r in{0,2,4,6}, tensor q copies of this eight-particle state and
a physical r-particle remainder. Every hopping term crossing a block boundary
changes the two block particle numbers, so its expectation is zero. Therefore

\[
E_0(N)\le q u_8+u_r.
\]

The remainder energies are recomputed directly from uniform U4,t1 CAR words and
bounded integer trial amplitudes. No submitted energy endpoint is accepted as
an input. The global state is specified by its block recipe and is not expanded.

At N=1,000,000, the new exact interval is

\[
-0.611636\le E_0/N\le u_8/8\approx-0.529475742405660.
\]

The width is approximately0.082160257594340 per site, a14.60% reduction from
the previous width0.096208555314006. The displayed approximate upper should be
rounded outward to-0.5294757424 when used as a decimal bound. A fresh replay
checks the six-site local lower certificate and upper tilings for N=10,12,14,64,
and1,000,000, covering every nonzero remainder. Its short-chain lower endpoint
uses only the overlap construction; existing tiling lower certificates can be
stronger for those short chains.

```sh
python -S results/marginal_graded_hubbard8/discovery/tiled_eight_upper.py
```

The corresponding inputs, source hashes, and exact intervals are under
`results/marginal_graded_hubbard8/tiled_eight_upper/`.

## Remaining work

This gives a reusable finite separator construction: identify a locally
admissible state, form its overlapping projectors, bound their joint Gram,
and turn the resulting incompatibility into a positive-operator certificate.
The witness selection and certificate discovery are not yet an efficient
general algorithm. A prescribed numerical probe of three alternative projector
witnesses and four penalties per witness gave no improvement over the accepted
fixed-witness certificate. The method still needs stronger physical upper states
with correlations across the current block boundaries,
tests on different interactions and geometries, and an accuracy-versus-cost
result as the allowed certificate complexity grows. None of these finite
results establishes general N-representability or a universal chemical compiler.

The exact replay driver is

```sh
python -S results/marginal_graded_hubbard8/discovery/projector_extendibility.py
```

Its inputs, fresh receipts, and source hashes are under
`results/marginal_graded_hubbard8/projector_extendibility/`. Focused and full
regression logs are recorded in the parent directory. Final validation passed
all432 regression tests in297.639seconds and all28 focused tests in7.166seconds.
An earlier425-test run caught an incorrect off-diagonal term in an intermediate
independent Gram test; that test was corrected to B n I−G and the full suite
rerun. Both the failed diagnostic and final passing logs are retained. The two
final standard-library drivers also replayed their exact certificates, and all
recorded source and input hashes matched.
