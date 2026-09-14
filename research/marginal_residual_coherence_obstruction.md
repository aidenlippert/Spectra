# Residual coherence directions and exact frozen-recipe limits

Two remaining positive-projector overlap witnesses yield independent stationary constraints beyond ENERGYv19/FAMILYv15. Exact three-source local mixtures limit the benefit of changing only these two coefficients to less than 2e-7 per site for W=0 and 8e-8 for W=1. These are bounds on a frozen-recipe subproblem; the full family gaps remain open.

## Independent constraints

The sources are |358>+|409> and |103>+|358>. Symmetry averaging and taking the reflection-odd part produce Y; T is the difference of its left and right embeddings. Removing the diagonal commutes with these operations. The full spin-word closure of the selected mixtures makes the removed moment exactly zero. Direct physical contractions reproduce the independently accepted full-density-matrix mismatches.

| Source case | Five-site entries | Six-site entries | Changed fermion bits | Overlap mismatch | Bound on T and open boundary sum |
|---|---:|---:|---:|---:|---:|
| W_zero | 8 | 56 | [8] | 0.000162442864213 | 1/2 |
| W_plus_1 | 16 | 128 | [2] | -0.000137849657358 | 1/4 |

The first source changes eight fermionic occupation bits. The second is a same-spin hop across four sites whose amplitude depends on the occupations of the other modes. It changes only two bits; a generic argument excluding all one-body support by bit count would be invalid. The physical Hamiltonian and its variable hopping profiles contain nearest-neighbor hopping, which these selected long-range entries exclude.

Exact row elimination finds particularly compact separating functionals:

- L0(M) = -16 M[413,1433].
- L1(M) = 8 M[1382,1433].

Both annihilate each of the 75 preceding offdiagonal operators exactly, including the two pure-coherence terms integrated in v19. Offdiagonal entries exclude all diagonal terms. Both lie in spin sector (4,2), outside the actual fixed half and charged projector supports, and neither is a nearest-neighbor transition. Their matrix on the two new operators is [[0,1],[1,0]], with determinant -1. Thus neither constraint is redundant modulo the full preceding affine family.

All Hermiticity, spin conservation, required signed symmetries, fermionic embedding translations and the six-site cyclic cancellation are checked exactly. Maximum absolute row sums bound Y by 1/4 and 1/8 respectively; the corresponding telescopes and open translated sums are bounded by 1/2 and 1/4. A positive-projector overlap mismatch violates stationarity of the particular local mixture. It does not violate that mixture’s local positivity or refute every mixture at its energy.

## Certified limits with the previous recipe frozen

For the complete fixed local matrix A, a positive trace-one mixture rho satisfying Tr(rho T0)=Tr(rho T1)=0 proves

`lambda_min(A + gamma0 T0 + gamma1 T1) <= Tr(rho A)`

for arbitrary real coefficients gamma0 and gamma1. Every older coefficient, physical profile and projector penalty remains fixed. Subtract the same penalty offset and divide by five. No fidelity inequality is required because the penalties are fixed. The existing accepted recipe at gamma=0 supplies the lower endpoint of this subproblem’s rigorous bracket.

| Case | Accepted periodic lower/site | Exact frozen-recipe ceiling | Maximum possible gain/site | Positive sources |
|---|---:|---:|---:|---:|
| W_zero | -0.64290658388107558 | -0.64290638460486316 | 1.99276212424e-07 | 3 |
| W_plus_1 | -0.66053609624452647 | -0.6605360247761094 | 7.14684171492e-08 | 3 |

The three-source witnesses have exact positive rational weights, exact trace one and both new moments zero. Their energies are recomputed from physical CAR actions, all old diagonal and offdiagonal corrections, and both fixed projectors using standard-library rational arithmetic. This includes the preceding pure-coherence coefficients; a separate nonzero offdiagonal test checks their contribution.

These ceilings do not bound joint reoptimization of older fields. The full ENERGYv19/FAMILYv15 gaps remain 4.49049126240775e-5 and 1.6702121170760045e-4 per site. Independence and a violated overlap alone therefore do not justify expecting useful improvement from a search that freezes the old recipe.

## Numerical probes and validation

| Case | Proposed new coefficients | Unaccepted proposed gain/site | Spectral evaluations |
|---|---|---:|---:|
| W_zero | ['7/1000000', '1/8000'] | 1.6e-07 | 9 |
| W_plus_1 | ['63/1000000', '-3/500000'] | 4e-08 | 14 |

These numerical energy proposals have no new production schema or exact PSD acceptance. The spectral probes check the active derivative and fresh physical matrix reconstruction and are capped at 250 optimization evaluations. Separate three-row LP searches use at most 20 pricing rounds and exact basis reconstruction. Independent physical replay, not floating-point convergence, establishes the ceilings.

All 40 focused tests passed, including 22 new tests and the 18 preceding coherence/cap tests. Coverage includes independent partial traces, exact diagonal-removal accounting, norm checks, dependence and injected-old-term refusal, physical nearest-neighbor exclusion, frozen old-coherence contributions, amplitude rescaling, malformed mixtures, nonzero moments and wrong version/scope refusals. One new test initially supplied profiles with invalid declared means and reflection; the existing validation correctly refused them. The fixture was corrected to valid nonuniform profiles and the failed version/log retained. No production source changed; the preceding full result remains 1128 tests plus 102 subtests and was not rerun for this diagnostic-only change.

The collector checks five current accepted receipts, 355 proof source-hash entries, and all 607 prior provenance files unchanged. It binds the current operators and ceilings to the accepted v19 energy and v15 mixture sources. Numerical proposals remain explicitly nonaccepting. No agents, GPU or paid resources were used.

No energy certificate or held-out transfer result changes in this turn. The preceding whole frozen recipe remains worse at its held-out target despite helpful individual terms. Joint reoptimization, other coherence directions, full-family attainment, general quantum representability, broader transfer and requested-accuracy scalability remain unresolved. The goal stays active.
