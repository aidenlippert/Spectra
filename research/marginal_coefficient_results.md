# Direct discovery in CAR coefficient space

The certificate search now operates on operator coefficients, without constructing a Fock basis or fixed-particle-number matrix. Five new certificates were discovered and replayed with exact rational arithmetic. The coupled eight-mode model admits a certified ground-energy interval of width 0.0000354093023573.

## Representation and mathematical contract

We represent a supporting energy bound for physical marginals by a positive operator polynomial modulo the fixed-number relation:

\[
H-bI=\sum_\alpha w_\alpha^\dagger Q_\alpha w_\alpha
+(\hat N-N)X+R,\qquad Q_\alpha\succeq0.
\]

The words obey the canonical anticommutation relations (CAR). Each Gram dictionary has a single particle-number charge. The multiplier X is a real Hermitian, number-conserving polynomial of ladder degree at most four. Its number-ideal term therefore vanishes on the target sector. This is a concrete, truncated dual representation; it is not a characterization of the entire physical-marginal cone.

Discovery matches independent real Hermitian coefficients through ladder degree six: 352 equations for six modes, 2,039 for eight. Numerical semidefinite optimization proposes b, X, and the Gram matrices. Export replaces the latter by rounded rational square factors. The checker reconstructs every coefficient, including Hermitian partners omitted from the numerical equation system, and charges the entire remaining residual:

\[
\eta=\sum_w |R_w|,\qquad E_0\ge b-\eta.
\]

This follows because every ladder monomial has operator norm at most one. Positivity follows from the exported square factors, independently of numerical solver status. No residual coefficients are discarded.

## Model and measured bounds

For M=2m modes at N=m particles, the matched model is

\[
H=\sum_{i<j\in A}n_i n_j+\sum_{i<j\in B}n_i n_j
-t\sum_{i=0}^{m-1}(a_i^\dagger a_{i+m}+a_{i+m}^\dagger a_i),
\]

where A and B are the two halves and t=1/5. Interaction strength sets the energy unit. This is an interacting finite fermion benchmark, not a molecular calculation.

| Modes | Dictionary | Certified lower bound | Variational upper, approximately | Interval width, approximately |
|---|---|---:|---:|---:|
| 6 | Full mixed cubic | 0.55099155545121 | 0.55100040032032 | 0.00000884486911 |
| 8 | Quadratic | 1.61050254103359 | 1.63786139471305 | 0.02735885367946 |
| 8 | Local cubic | 1.61050092122707 | 1.63786139471305 | 0.02736047348598 |
| 8 | Full mixed cubic | 1.63782598541069 | 1.63786139471305 | 0.00003540930236 |

Upper bounds and widths in the table are rounded displays; exact fractions are saved in the certificates and replay output. The eight-mode mixed-cubic upper bound is exactly 8189306801722658/4999999895081123.

The asymmetric six-mode control changes matched hopping weights to 1, 7/10, 13/10 and adds an edge (0,4) of weight 2/5. Its newly certified lower bound is 0.54086857742176. This run does not supply a new asymmetric upper witness.

The quadratic dictionaries contain linear, pair, and particle-hole words. The local-cubic extension adds within-half triples and within-half cubic words mixed with linears. The full mixed extension allows all words of the form a_k† a_j a_i in its charge-minus-one block and their adjoints in the opposite block. Its pure triple blocks remain local. All Gram entries within each block are independently optimized.

The full mixed dictionary yields an approximately 773-fold narrower eight-mode interval than the quadratic dictionary. Local cubic terms do not improve the achieved certified bound in this run; the tiny decrease is within the export/solver accuracy difference. This comparison does not establish which individual cross-group words are necessary.

## Upper witness without configuration enumeration

A valid variational state has one particle per matched pair, with integer amplitude a_k shared by every flavor-ordered configuration containing k particles on the left. Its norm and energy are evaluated by binomial counts:

\[
Z=\sum_k\binom mk a_k^2,
\]
\[
U=\frac{\sum_k\binom mk[\binom k2+\binom{m-k}2]a_k^2
-2t\sum_{k=0}^{m-1}\binom mk(m-k)a_k a_{k+1}}{Z}.
\]

Discovery chooses amplitudes through a small symmetric variational matrix of dimension m+1 and rounds them to integers. Exact replay recomputes the quotient and verifies that the certificate Hamiltonian equals the matched model. This is an upper bound to the full N-particle ground energy; it does not assume that the ansatz contains the ground state. Independent tests compare the combinatorial formula with signed fermionic action for four, six, and eight modes, including negative hopping and alternating amplitude signs.

## What the computation establishes

Both discovery and certificate verification now avoid explicit occupation-sector matrices. Numerical eigendecompositions remain for Gram factorization and the small variational ansatz. Discovery still requires NumPy, SciPy, and CVXPY; exact interval replay uses only the Python standard library.

The eight-mode mixed search uses Gram dimensions [8,8,28,28,64,8,8,232,232], 443 multiplier basis elements, and 50,004 nonzero entries across its Gram-to-coefficient maps. Its exported certificate has 480 factor rows and 32,556 nonzero factor entries. The recorded build and solve times were 2.47 and 14.17 seconds respectively in one run; these are observations, not a scaling or speedup result.

The eight-mode residual allowance is 0.00003539152831, accounting for almost all the certified interval width. The rational proposed b lies about 1.78e-8 below the upper witness. Thus this run separates two next questions:

1. Can higher-accuracy solving and rational reconstruction reduce the residual allowance while preserving the dictionary?
2. Can a smaller, transferable operator dictionary retain the full mixed bound? The largest dense blocks have dimension M+M*binom(M,2), so removing Fock enumeration alone does not make this approach cheap.

No result here establishes a universal compact representation, efficient general N-representability, or general chemistry prediction. The concrete advance is direct discovery and exact verification of tight dual certificates on interacting six- and eight-mode examples.

## Reproduction and validation

From the workspace root, using the numerical Python environment:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_coefficient --modes 8 --t 1/5 --family mixed
```

Independent exact replay, without site packages:

```sh
python3 -S -m experiments.marginal_collective results/marginal_coefficient/m8_t1_5_mixed.json
```

All 27 tests across test_marginal_hunt_car, test_marginal_hopping, test_marginal_symbolic, test_marginal_coefficient, and test_marginal_collective passed. All five new lower certificates and four matched upper witnesses replayed successfully in a `python3 -S` process that checked NumPy, SciPy, and CVXPY were not loaded.

Implementation: `experiments/marginal_coefficient.py` and `experiments/marginal_collective.py`. Exact results and source/artifact hashes: `results/marginal_coefficient/replay.json` and `results/marginal_coefficient/manifest.json`. Existing historical experiment outputs were preserved.
