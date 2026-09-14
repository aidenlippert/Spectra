# Physical-marginal hunt: findings, exact witness, and next attack

Four independent agents investigated certificate algebra, marginal geometry,
prior art, and obstructions. The root reviewed and corrected their notes and
implemented an exact separation example. This is a mathematical research result
and an executable algebra diagnostic, not a trained compiler or molecular solver.

## The representation to pursue

Learn a small, orbital-adapted dictionary of fermionic operators that nearly
annihilate a physical trial state **and can jointly reproduce the Hamiltonian
as a positive Gram form**. Its dual image supplies the supporting inequalities
of the physical-marginal body.

For a dictionary column w, search

    H - b I = w† Q w + Z + R,
    Q >= 0,  P_N Z P_N = 0,  ||P_N R P_N|| <= eta.

P_N projects onto the declared particle-number sector. All operators and
residuals in this equation are Hermitian where required. The lower bound is
b-eta. A separately certified variational energy u gives [b-eta,u]. The Gram
factor Q=L†L supplies the actual squares B_alpha=(Lw)_alpha.

At zero residual, for a normalized physical state Psi,

    E_Psi - b = sum_alpha ||B_alpha Psi||².

Thus tight certificates require near-annihilators. With residual R the certified
width is sum ||B_alpha Psi||² + <R> + eta, plus any upper-evaluation error.
Finding a kernel alone is insufficient: the same operators must match H.
This complementary-slackness mechanism is also explicit in Hastings' SOS
analogue of Wigner's rule, equations (2)-(3) of his
[2024 paper](https://arxiv.org/html/2412.03564v1).

For a fixed dictionary coefficient matching is linear in Q and b; dictionary
selection is the nonlinear learning problem. Impose homogeneous charge blocks.
Cubic-word Gram entries generally require a 3-RDM or direct state contractions.
The 2-RDM alone does not supply them.

**Open conjecture:** chemically useful directions admit transferable small
dictionaries whose certified energy gap remains small through strong correlation.
This does not follow from convexity or from the existence of separators.

## Exact algebraic mechanism

If B is a homogeneous-charge odd polynomial of ladder degree at most three,

    C(B) = B†B + BB† >= 0

has ladder degree at most four after CAR normal ordering. Its degree-six symbols
cancel because two odd cubic symbols anticommute. It therefore defines a valid
2-RDM inequality. Pure triple annihilators give T1-type conditions; mixed cubics
give T2-type conditions. A single cubic square generally retains degree six.

This is existing positivity algebra, not a newly discovered family. Its practical
use here is selective generation: propose a small coefficient vector B for which
the current pseudo-marginal has negative expectation of C(B), verify, and add it.
Once full DQG holds, more quadratic-square cuts cannot strengthen it. Likewise,
once the full relevant T family holds, selecting more members cannot improve it.

Pairing is not automatically tight: both B Psi and B† Psi contribute. For a bare
canonical mode their squared norms sum to one. General dressed SOS constructions
are needed beyond merely selecting anticommutator cuts.

## An exact scalable example beyond DQG

Take 2r disjoint triples T of spinless modes, with 6r modes and N=3r particles.
Use unit interaction coefficients and

    H_r = sum_T sum_(i<j in T) n_i n_j,    S_T = sum_(i in T) n_i.

This diagonal density-interaction model is deliberately analytically soluble.
It diagnoses the strength of a certificate family without a solver or learning
confound. It is not a realistic chemical benchmark.

Define a diagonal pseudo-marginal by

    <n_i> = 1/2,
    <n_i n_j> = 1/8 inside a triple, 1/4 across triples,

with off-diagonal pair-basis 2-RDM entries zero. Antisymmetric permutations are
understood; the canonical unordered-pair basis avoids double counting.
The fixed-N contraction holds for every i:

    sum_(j != i) <n_i n_j>
       = 2/8 + (6r-3)/4 = (N-1)/2.

### Full DQG feasibility

D is the Gram matrix of pair annihilators a_j a_i, i<j. It is diagonal with
entries q_ij in {1/8,1/4}. Q is the Gram matrix of their adjoints, diagonal with
entries 1-p_i-p_j+q_ij=q_ij. Both are PSD.

G is the Gram matrix of a_i† a_j, with entries <(a_i†a_j)†(a_k†a_l)>.
It splits into:

- A number-operator block P with diagonal 1/2, off-diagonal entries 1/8 within
  a triple and 1/4 across triples.
- Independent transition entries p_j-q_ij, equal to 3/8 or 1/4.
- Zero cross blocks and zero off-diagonal transition entries.

The number block has the exact factorization

    P = J_(6r)/4 + blockdiag_T[(3 I_3 - J_3)/8].

Each 3 I_3-J_3 is a graph Laplacian, the sum of the three outer products
(e_i-e_j)(e_i-e_j)^T. Hence P>=0. The augmented identity/number moment matrix
is also PSD because P-pp^T is precisely that block-diagonal Laplacian.
The one-body and hole matrices are I/2, and the number variance is zero.
The full G constraints, rather than just its diagonal, therefore hold.

### Exact gap and its closure

The candidate has energy 3r/4. Conversely, the G-compatible square inequality

    sum_T <(S_T-3/2)²> = 2<H_r> - 3r/2 >= 0

proves that the DQG optimum is exactly 3r/4. The constant in each square is
allowed because I=Nhat/N on this sector (or via the augmented moment matrix).

For B_T=a_i a_j a_k on a triple, exact CAR gives

    {B_T,B_T†} = 1-S_T + sum_(i<j in T) n_i n_j.

Its pseudo-expectation is -1/8. Summing these physical inequalities gives

    sum_T {B_T,B_T†} = H_r-r I    on N=3r.

Therefore H_r>=r. Occupying two modes in r triples and one mode in the other r
triples gives a physical product state with energy r. The exact ground energy is
r, and 2r anticommutator cuts (4r squares) close the entire r/4 DQG gap.

**Result:** a linear number of constant-support cubic certificates removes an
extensive error that survives the full DQG relaxation. The construction is an
elementary occupancy/T1 example; no novelty or quantum-chemistry claim is made.

### Geometric extension and the next obstruction

Any unitary orbital rotation d=Ua preserves CAR and particle number. Replacing
a by d in H_r and every B_T preserves the identity and the bounds. A dense
orbital basis can make the coefficient arrays dense while the certificate
remains a product of three linear orbital forms. The common U costs O(M²) to
store. Expanding those forms and contracting a generic Hamiltonian still costs
work: compact factorization is not a claim of free evaluation.

This illustrates why learning the orbital representation could matter for a
sparse dictionary. It does not strengthen the full basis-invariant SOS cone,
and the rotated family remains unitarily equivalent to an easy diagonal model.

Adding hopping between triples is the concrete next obstruction: local occupancy
certificates will generally cease to be tight. A norm-bounded residual keeps the
old proof sound, but its error may be extensive. The useful question is whether
adaptive mixed cubic/linear words close this new gap with a small dictionary.

## What the source hunt changed

- [Mazziotti 2021](https://arxiv.org/abs/2103.17155) already formulates a
  dual-cone 2-RDM calculation. The cone language is established.
- [Hastings 2022](https://arxiv.org/abs/2205.12325) identifies degree-six
  fragments that recover perturbative information missed by degree four.
- [Hastings 2024](https://arxiv.org/html/2412.03564v1) constructs self-consistent
  SOS witnesses with favorable model results, while identifying strong diagonal
  chemical interactions as a limitation of that implementation.
- [Hastings June 2026](https://arxiv.org/html/2606.31765v1) addresses selected
  large density-density terms and related interactions. This is a concrete
  algebraic starting point, not evidence of a solved universal compiler.
- [Delgado-Granados and Mazziotti 2026](https://arxiv.org/html/2603.05524v1)
  train an ICNN penalty for marginal geometry inside DQG optimization. Their test
  split includes geometries of the test molecule in training. A learned penalty
  by itself does not certify a supporting inequality or a ground-energy bracket.

The strongest attack emerging from this pass is **learned dressed operator
dictionaries with exact Hamiltonian matching**, seeded by selective T constraints
and the physical state's approximate kernels. Compare against full available
T conditions, deterministic cut selection, active-space constraints, and the
self-consistent SOS construction before attributing any gain to learning.

## Verification receipt and limits

Executed using Python 3.13 standard library only:

    python3 -m unittest tests.test_marginal_hunt_car -v
    python3 -m experiments.marginal_hunt_car
    python3 -m experiments.marginal_hunt_witness

Six focused tests pass. They include independent occupation-basis evaluation of
all 1,296 length-four words on three modes (10,368 state-action comparisons),
pure/mixed cubic and cubic-plus-linear anticommutators, a deliberately corrupted
identity, all 1,746 D/Q/G entries for six modes, and exact fixed-N arithmetic.
The six-mode diagonal ground energy is independently enumerated. Larger receipt
rows through 768 modes instantiate the analytic formulas; they are **not**
768-mode SDP solves, dense diagonalizations, or measured scaling benchmarks.

The CAR code handles real rational coefficients in this diagnostic; it is not
a general complex-coefficient or interval-certified production verifier. No
Lean formalization, chemical integrals, molecular energy calculations, learning,
or runtime advantage have been established. Historical experiments were not run
or modified; the new code uses no existing scientific dependencies.

Artifacts: [CAR code](../experiments/marginal_hunt_car.py),
[witness replay](../experiments/marginal_hunt_witness.py),
[tests](../tests/test_marginal_hunt_car.py),
[CAR results](../results/marginal_hunt_car.json),
[witness results](../results/marginal_hunt_witness.json).

Supporting notes: [algebra](marginal_hunt_algebra.md),
[geometry](marginal_hunt_geometry.md), [prior art](marginal_hunt_prior_art.md),
[obstructions](marginal_hunt_obstructions.md). Where earlier exploratory
phrasing is less precise, this synthesis states the reviewed conclusions.
