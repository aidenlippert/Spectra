# Compact coherence obstructions and a fixed-recipe limitation

Two exact, independent stationary constraints remain after full diagonal spin-word closure. They can be made purely offdiagonal and sparse. However, exact local-mixture ceilings show that changing only these two coefficients, while freezing every older coefficient, cannot materially improve the current energy recipes.

## Exact new constraints

Start from the normalized positive projectors onto |346>+|409> and |314>+|614>. Average under particle-hole and spin flip, take the reflection-odd part Y, and form T=Y_left−Y_right. Remove Y’s diagonal: the complete spin-word constraints make its expectation zero in the selected symmetry-averaged mixtures. Direct contractions of the resulting pure-coherence T still equal the independently reconstructed positive-projector overlap mismatch.

| Source case | Five-site entries | Six-site entries | Changed fermion bits | Exact mismatch, decimal approximation | Norm bound on T |
|---|---:|---:|---:|---:|---:|
| W_zero | 8 | 60 | [4] | -0.000219649124863 | 1/2 |
| W_plus_1 | 16 | 120 | [6] | -0.000270921208986 | 1/2 |

Each Y has maximum absolute row sum 1/4, giving spectral norm at most 1/4. Both T and any open translated telescoping sum therefore have norm at most 1/2. Exact signed-permutation checks establish Hermiticity, spin conservation, all required symmetries, fermionic left/right translation and zero periodic sum. The original witnesses are positive projectors; their nonzero overlap difference is a consistency violation, not negative local positivity. The pure-coherence differences themselves need not be positive.

The four-bit constraint is an occupation-dependent spin exchange and shares support with an old spin term. Merely counting changed modes would not establish independence. The exact linear functionals used instead are:

- L0(M)=M[1370,1433]−M[350,413].
- L1(M)=M[1337,1637].

Both functionals annihilate all 73 existing nondiagonal operators. All selected entries lie in spin sector (4,2), excluding the actual fixed half and charged projectors. Four or six changed bits exclude physical one-body hopping, and offdiagonality excludes every diagonal term, including all 120 spin-word directions. On the two new T operators, the functional matrix is diag(1/8,1/16), with determinant 1/128. This proves two independent directions modulo the entire ENERGYv18/FAMILYv14 affine family.

## Exact limitation when older coefficients are frozen

Let A be the complete accepted local matrix, including every frozen old coefficient and both fixed projector penalties. For a positive trace-one local mixture rho with Tr(rho T0)=Tr(rho T1)=0,

`lambda_min(A + gamma0 T0 + gamma1 T1) <= Tr(rho A)`

for every real pair gamma. Subtracting the unchanged projector penalty cost and dividing by five gives the following exact ceilings. The witnesses do not require projector fidelity inequalities because both penalties remain fixed. All source vectors, weights, moments and local energies are replayed with standard-library rational arithmetic and physical CAR actions.

| Case | Accepted periodic seed lower/site | Fixed-old-recipe ceiling | Maximum possible gain | Positive sources |
|---|---:|---:|---:|---:|---:|
| W_zero | -0.64291211911474899 | -0.6429120814230429 | 3.76917060319e-08 | 1 |
| W_plus_1 | -0.66056268336259027 | -0.66056251056563819 | 1.72796952167e-07 | 3 |

These bounds cover unrestricted real new coefficients. They prove a narrow limitation of the fixed old recipes, not of joint reoptimization. The full previous-family gaps remain about 5.75e-5 and 3.80e-4 per site; neither full-family numerical limit is resolved here.

The bounded two-variable spectral probes found zero coefficients and no gain for W=0. W=1 proposed coefficients −93/500000 and 179/500000, with a proposed gain of 1.4e-7/site. This remains an unaccepted numerical energy proposal: no new energy schema or exact PSD lower replay was introduced. Both probes check an active spectral derivative and fresh physical matrix reconstruction. The exact ceilings above independently bound what any two-coordinate search could achieve.

A separate three-row LP search proposed the positive mixtures. It used at most 20 pricing rounds, bounded integer amplitudes and an exact fraction-free reconstruction before the independent physical replay. Numerical success is not acceptance. The initial replay implementation misread the existing action dictionary as a list and failed; focused tests caught the same error. That version and its failed logs are preserved. The corrected dictionary indexing passed the independent occupation-energy tests and both exact replays.

## Validation and remaining scope

All 39 focused tests passed: 30 overlap/coherence tests and 9 fixed-recipe cap tests. They include independent partial traces, norm checks, diagonal-removal accounting, amplitude rescaling, direct vacuum occupation energies, nonzero-moment refusal, malformed mixtures and scope refusals. The collector audits five accepted receipts and 350 proof source-hash entries; all 548 prior provenance files are unchanged. No production source changed and no full suite was rerun. The preceding full result remains 1074 tests plus 102 subtests.

No energy lower bound or held-out transfer result changes in this turn. The new constraints are ready for a joint reoptimization experiment, where older diagonal, offdiagonal, profile and penalty coefficients can respond. Other coherence directions may also matter. General quantum representability, full-family numerical attainment, broader molecular/long-range/higher-dimensional transfer, and scalability at requested accuracy remain unproved. No agents, GPU or paid resources were used; the goal stays active.
