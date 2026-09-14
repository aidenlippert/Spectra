# Hopping attack: adaptive operators improve bounds, paired squares hit a limit

## Result

The original diagonal model now has inter-triple hopping. Deterministic selection
of mixed cubic operator directions improves certified ground-energy intervals in
all four tested cases. Random directions at the same eight-pair budget barely
improve the original local certificates. A larger dictionary with independently
weighted charge blocks produces much tighter certificates than a dictionary
restricted to equal weights for conjugate pairs.

This is a finite-sector experiment in six modes and three particles, not a
trained model, a molecular calculation, or evidence of scalable quantum solving.
The exact verifier uses the 20-dimensional occupation sector. All 34 exported
certificates replay, and 13 focused tests pass (six prior CAR tests, seven new
hopping tests including a solver/export integration check).

## Hamiltonians

The interaction coefficient is one, defining abstract energy units:

    H(t) = sum_(i<j in A or B) n_i n_j
           - t sum_(i,j) w_ij (a_i† a_j + a_j† a_i),
    A={0,1,2}, B={3,4,5}, N=3.

Matched edges are (0,3),(1,4),(2,5), with all weights one. The asymmetric case
uses weights 1,7/10,13/10 on those edges and an additional (0,4) edge of weight
2/5. This changes the hopping symmetry but is not an unseen chemistry family.
Both are tested at t=1/5 and t=1. These values and the single random seed were
chosen for an exploratory diagnostic, not a preregistered statistical study.

## What was actually learned or searched

No neural network was trained. The adaptive method reads the numerical dual
matrix of the current SOS coefficient-matching problem. It constructs the moment
matrix for a candidate operator dictionary and selects its most negative
anticommutator direction. One direction and its adjoint are added per round,
for eight rounds. The Hamiltonian is then matched again with PSD Gram matrices.

The mixed dictionary contains 96 words per charge: six linear words and 90
words a_k† a_j a_i (i<j), together with the adjoint dictionary. Density-assisted
words n_j a_i form a smaller 36-word dictionary. Full pure triple dictionaries
are excluded: with N=3 those alone could encode arbitrary sector matrices.

Each selected direction is a generally dense coefficient vector over the pool.
Eight directions means eight linear combinations, not eight monomials. The
baseline quadratic Gram blocks remain present. No asymptotic memory or runtime
advantage is claimed.

The physical upper witness is a rounded numerical exact-diagonalization vector.
Its rational Rayleigh quotient is evaluated exactly. It is not used to select
the adaptive directions. This experiment tests the dual search; it does not yet
implement primal-guided near-annihilator learning or transfer across systems.

## Certified widths at equal eight-pair budget

Widths are upper minus verified lower, including factor-rounding and numerical
matching residuals. They are not the numerical SDP objective error.

| Hamiltonian | t | Quadratic + original local cubics | Random eight pairs | Selected eight pairs |
|---|---:|---:|---:|---:|
| Matched | 0.2 | 0.0409034 | 0.0409034 | 0.0258776 |
| Matched | 1.0 | 0.00330625 | 0.00330861 | 0.000102867 |
| Asymmetric | 0.2 | 0.0397689 | 0.0397817 | 0.0222958 |
| Asymmetric | 1.0 | 0.00315364 | 0.00316170 | 0.0000746152 |

The selected search improves the stronger-hopping widths by approximately
32-fold and 42-fold relative to the local-cubic baseline. These ratios concern
interval width, not runtime. One random seed supports only the reported
comparison, not statistical superiority over random selection generally.

The narrower density-assisted dictionary fails to close the gap: its complete
pool leaves widths approximately 0.04090 and 0.00330 in the matched cases. Its
initial unpaired direction-selection attempt also largely stalls. Those null
results are retained in summary.json rather than discarded.

## The algebraic distinction exposed by the control

For t=0.2, compare the full mixed dictionary in two forms:

    Paired:    w_-† Q w_- + w_+† Q w_+,       Q >= 0
    Untied:    w_-† Q_- w_- + w_+† Q_+ w_+,   Q_-,Q_+ >= 0.

Here w_+ consists of the adjoints of w_-. Real symmetric Gram matrices are used.
The paired construction is an anticommutator sum. Its degree-six cancellation
is automatic. Untied weights are more general; matching and any higher-order
residual must be checked on the declared sector.

| Hamiltonian, t=0.2 | Full paired dictionary width | Full untied dictionary width |
|---|---:|---:|
| Matched | 0.0213302 | 0.00000500379 |
| Asymmetric | 0.0198451 | 0.0000135652 |

This demonstrates substantially stronger **achieved certificates** from the
untied construction. Numerical optimization suggests a limitation of the paired
family on these cases; we have not independently certified an upper bound on
the best possible paired-SOS objective, so this is not a theorem of optimal
separation between the two proof systems.

The full untied lower/upper intervals are:

    matched:    [0.550995396532, 0.551000400320345...]
    asymmetric: [0.540863615947, 0.540877181133642...].

The factors are exact rational finite-sector witnesses. We did not establish
that their CAR expansion reduces to a compact two-body identity using a short
number-ideal certificate. That remaining step is essential for scaling beyond
the enumerated sector.

## Degenerate directions and subspace selection

At the matched local-cubic solution for t=0.2, the most negative candidate
anticommutator eigenvalue is about -0.0204165 with numerical multiplicity 12.
There are 18 negative eigenvalues in total. Selecting the whole negative
subspace in one pass produces width 0.0223794; a random subspace of equal
dimension leaves width 0.0409144. In the asymmetric case the corresponding
widths are 0.0216812 and 0.0397809.

Recycling such subspaces improves the matched width to 0.0203593 with 36
directions per charge. The asymmetric run reaches 0.0168667 at 32 directions,
then worsens to 0.0174906 after further enrichment because the numerical solve
and rational residual become less accurate. These are exploratory extra-budget
runs. They do not establish that selecting a subspace solves the compactness
problem, and their budgets differ from the eight-pair table.

## Exact verification, including failed proposals

The numerical optimizer only proposes Gram factors. For each PSD block we
factor its symmetrized numerical Gram matrix, clip negative eigenvalues, and
round the factor entries to integers divided by 10^6. Positivity is then
defined by the exact factor product L^T L, not assumed from solver status.

The verifier independently builds integer occupation-basis action matrices W_i,
constructs every rational square F_alpha^T F_alpha, and computes

    R = H - b I - sum_alpha F_alpha^T F_alpha,
    eta = max_i sum_j |R_ij|,
    lower = b-eta.

For a symmetric residual the row norm bounds the operator norm. The nonzero
integer trial vector x gives the exact rational upper bound x^T H x/(x^T x).
Float Hamiltonians, malformed factors, wrong sectors and invalid trial vectors
are rejected. A bad but well-formed scalar or factor cannot fool the check:
it increases the residual charge and weakens the safe bound.

CLARABEL failed on the first redundant sector-matching formulation; SCS supplied
the recorded proposals. Some enriched solves reported possible inaccuracy.
Neither event invalidates the separately replayed bounds, but it affects their
tightness. Several numerical estimates are therefore deliberately not described
as exact SDP optima. A sound factorization and a small solver residual are
different claims.

## Analytic reference for the matched model

The matched model conserves each flavor population n_i+n_(i+3). In the
(1,1,1) population block it reduces to three flavor spins. The spin-3/2 block
has eigenvalues

    2-t +/- sqrt(1+2t+4t²),
    2+t +/- sqrt(1-2t+4t²).

The spin-1/2 blocks and the (2,1,0) flavor-population permutations contribute
1-t and 1+t, each with combined multiplicity eight. This lists all 20 levels.
For t>=0 the minimum is

    E0(t) = 2-t-sqrt(1+2t+4t²).

For example it is below 1-t because sqrt(1+2t+4t²)>=1, and below the other
parity minimum because both 2t and the difference of the two square roots are
nonnegative. This checks the numerical reference without asserting a compact
SOS factorization: pulling arbitrary spectral factors back to fermionic
polynomials can require high-degree projectors.

## The next mathematical target

Search for **small independently weighted operator families whose combined
high-degree terms vanish modulo fixed-N identities**, using a compact symbolic
verifier rather than a full occupation-sector matrix. Automatic anticommutator
cancellation is a useful seed, but imposing it on every pair restricts the
search. The variables to learn are the operator subspaces and the relationships
among their Gram weights.

The immediate unsolved obligation is to extract that short symbolic identity
from the successful finite-sector factors. Until then the experiment establishes
stronger finite certificates and a concrete restriction to relax, not the
AlphaFold moment of quantum matter.

## Reproduction and artifacts

Use the existing environment at
/opt/homebrew/Caskroom/miniconda/base/bin/python with NumPy, SciPy and CVXPY.
No new package was installed or added to the historical project requirements.
The module records numerical proposals, while replay is exact rational.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_hopping --stage baseline
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_hopping --stage paired
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_hopping --stage random
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_hopping --stage full
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_hopping --stage anticommutator
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_hopping --stage subspace
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_hopping --stage recycled
    python -m experiments.marginal_hopping --stage replay
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest tests.test_marginal_hunt_car tests.test_marginal_hopping -v

The additional stages were initially run through equivalent direct module calls
and then given named CLI entry points. Re-running may change floating optimizer
details; saved rational certificates replay without re-solving. Some early control
filenames differ in word order from the consolidated runner's filenames.

[Implementation](../experiments/marginal_hopping.py),
[tests](../tests/test_marginal_hopping.py),
[34-certificate replay](../results/marginal_hopping/replay.json),
[runtime and hashes](../results/marginal_hopping/manifest.json),
[selected-direction results](../results/marginal_hopping/paired_summary.json),
[full mixed results](../results/marginal_hopping/full_mixed_summary.json),
[paired-family control](../results/marginal_hopping/anticommutator_summary.json).
