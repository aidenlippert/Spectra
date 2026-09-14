# Certified energy transfer for quadratic boundary filters

The quadratic-filter construction now has a physical energy certificate for
the million-site, half-filled open Hubbard chain at U=4,t=1. The exact
certificate gives the outward-rounded ground-energy interval

\[
-0.611636\le E_0/N\le-0.56726244870.
\]

Its width is approximately 0.0443735512909309 per site, a 16.96% reduction
from the preceding adapted linear-filter interval. The new upper improves
by approximately 0.00906227771539726 per site. The lower bound still comes
from the independently replayed six-site all-Fock certificate.

The accepted gate is

\[
F=I-a h+b h^2,\qquad a=0.356039,\quad b=0.177774,
\]

where h is the physical hopping on one interface. These are exact rational
parameters in the replay. The building block is the original validated H8
state; this result does not require the adapted source from the preceding
milestone. Every interface receives the same even, number-conserving gate.

The two endpoints of the transfer calculation enclose this trial state's
energy. They are not both bounds on the ground energy: only the trial
energy's upper endpoint supplies a variational ground-energy upper bound.

## Retaining the correlations that scalar merging lost

The earlier exact counterexample showed that quadratic filters change
remote endpoint observables. The new construction therefore contracts the
actual boundary tensor instead of assuming those observables are inherited.

For a real block amplitude phi(l,m,r), with single-site occupations l,r
and middle configuration m, define

\[
G_{(l,l'),(r,r')}=\sum_m\phi(l,m,r)\phi(l',m,r'),
\]

\[
B_{(r,r'),(l,l')}=\langle r',l'|F^\dagger F|r,l\rangle.
\]

Each paired single-site index has dimension 16. With
e_(l,l')=delta_(l,l'), the exact q-block norm is

\[
\mathcal N_q=e^T G(BG)^{q-1}e.
\]

The implementation computes G from an explicit physical integer block
state, rather than accepting an arbitrary supplied marginal. Occupation
indices stay in physical site order throughout. For adjacent hopping,
the preceding Jordan-Wigner strings cancel as an operator identity;
no state-dependent block-parity correction is omitted.

## The energy terms that touch a filter

For blocks of length L>=4, partition the Hamiltonian into:

- A middle operator for each block: onsite terms on sites 1 through L-2
  and hopping bonds 1 through L-3, using zero-based labels.
- The left and right outer endpoint onsite terms and adjacent hoppings.
- One four-site patch per interface: the onsite terms on its two contact
  sites, the hopping across the interface, and the internal hopping next
  to each contact.

This assigns every term exactly once. The middle operator is disjoint from
all contact filters. Each four-site patch touches only its own filter.
Its contribution is therefore computed using F-dagger H_patch F, with all
other filters represented by their F-dagger F kernels.

Let J,L,R denote block tensors with middle, left outer, and right outer
insertions. Let W denote a two-block tensor with the dressed four-site
patch inserted and the two outer endpoints retained. W is computed from
physical three-site partial traces of each block. It is a 16-by-16 tensor;
the code never constructs a sixteen-site wavefunction.

Write T=BG. A row carrying the norm A_n, its preceding value P_n=A_(n-1),
and the accumulated energy E_n evolves by

\[
A'=AT,\qquad P'=A,\qquad
E'=ET+A(BJ)+P(BW).
\]

The two-block seed explicitly includes the first middle and left-boundary
insertions, the second middle insertion, and the first dressed patch. The
final contraction adds the right boundary through P_n B R e. A separate
one-block branch directly evaluates its full Hamiltonian.

This gives a 48-dimensional augmented linear transfer. Its physical input
is a bounded block state and a local gate, not a global many-body vector.

## Rounding without expanding the normalization

Exact rational recurrence is available for up to 64 blocks. Longer chains
use interval matrix powers with fixed precision:

1. Each rational entry is enclosed on a dyadic grid.
2. Each product uses all four interval-endpoint products, followed by
   outward floor and ceiling operations.
3. Matrices are rescaled by a common positive power of two. That scale is
   tracked symbolically through every multiplication.
4. The final norm must have a strictly positive lower endpoint. Otherwise
   the verifier refuses the energy certificate.
5. The shared scale cancels in the energy-to-norm ratio. Taking the minimum
   and maximum of its four endpoint quotients encloses the trial energy.

At 125,000 blocks, the 160-bit calculation uses 16 matrix squarings and
8 row-times-power products. Its energy-density enclosure width is
approximately 9.84883684481505e-38. The final common binary exponent is
82,969; it is stored as an integer, not expanded into an enormous norm.

This numerical enclosure width is distinct from the much larger
ground-energy interval width. It establishes that rounding is negligible
for this certificate, not that the trial state is nearly exact.

The transfer size is fixed for this chosen block and gate family. The
explicit block construction still uses a large part of the finite H8
sector. The method does not establish efficient arbitrary block discovery,
general N-representability, or a precision-versus-cost guarantee for all
Hamiltonians. Refusal at inadequate precision remains possible.

## Independent checks and reproduction

The replay performs four different checks:

- Direct physical CAR calculations verify norms and full Hubbard energies
  on one, two, and three four-site blocks, including a source with charge
  fluctuations and mixed linear/quadratic gates.
- Exact finite transfer energies lie inside the rounded enclosures.
- At one million sites, the linear-filter specialization encloses the
  exact energy from the earlier, independently implemented scalar merging
  formula.
- A fresh six-site PSD proof supplies the ground-energy lower bound, and
  a fresh adapted-block replay supplies the previous upper comparator.

Signed interval-product tests include both very large and very small
scales. A physical product annihilated by I-h^2 is rejected by both the
exact and rounded contractions. The 96-bit and 160-bit large-chain
enclosures also overlap.

```sh
OPENBLAS_NUM_THREADS=1 python -S results/marginal_graded_hubbard8/discovery/boundary_transfer.py
```

The fresh replay completed in 53.147 seconds, including the lower PSD
proof. Inputs and exact receipts are in
`results/marginal_graded_hubbard8/boundary_transfer/`. Numerical discovery
is a separate script; its output is never trusted by the accepting replay.
The 32 focused tests passed in 34.103 seconds. All 457 regression tests
passed in 321.626 seconds. All 14 recorded source/input hashes matched
afterward. Validation is recorded in `results/marginal_final_validation.json`.

## Remaining work

The scalar-closure obstruction has been bypassed with a verified energy
contraction that retains boundary correlations. The accepted a,b values
were numerically proposed; their global optimality is unproved. The block
state and broader contact-operator families can now be optimized against
the full transfer energy. Stronger transferable lower certificates remain
necessary to narrow the remaining ground-energy interval. The current
transfer compiler handles the specified Hubbard model; molecular Coulomb,
thermal, dynamic and synthesis claims remain unsupported.
