# Spectra tensor response constructor

Spectra can now construct a correlated response directly in tensor coordinates,
grow those coordinates using the physical residual, and independently check the
result without enumerating the many-body sector. This is a working numerical
construction backend. It does not require an LLM, a pretrained network, a full
wavefunction teacher, or the previous determinant/Krylov pool.

## Run

Use Python 3.10 or later, NumPy, SciPy, a C11 compiler, and GMP development headers
and library. The accepting process itself imports no numerical Python packages.

```sh
python -m pip install -r spectra_tensor/requirements.txt
python -m spectra_tensor request --sites 20 --omega 2 -1 8 --out request20.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m spectra_tensor solve \
  --request request20.json --out tensor20
```

Output directories must be new. Every checked attempt retains its candidate,
request, exact receipt, solver measurements, and acceptance decision. Exit code
2 means the requested accuracy was not achieved inside the supplied bond and
sweep budget. It is not reported as success.

A pure-Python integer backend is available with `--backend python`; use it for
small cross-checks. The GMP backend compiles the included C arithmetic kernel
into a private temporary cache. `CC`, `CPPFLAGS`, and `LDFLAGS` are passed as
argument lists, without a shell. An installation outside standard include and
library paths can be selected with `SPECTRA_GMP_PREFIX=/path/to/gmp/prefix`.
Linux execution was tested in this delivery; macOS execution was not available.

Replay any individual candidate in a separate process:

```sh
python -B -S -m spectra_tensor verify \
  tensor20/query_00/bond_096.candidate.json \
  --request tensor20/query_00/request.json --out replay20.json
```

Run the package tests from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m unittest \
  discover -s tests -p test_tensor_response.py -v
```

## Physical statement

The supported model is a finite, repulsive, spin-independent Hubbard Hamiltonian
with real rational hopping coefficients and onsite interactions. The source is
an alternating-spin product state at even half filling. All results concern

`g(z) = b^T (z I - H)^(-1) b`, with `z = omega + i eta`, `eta > 0`.

This is a scalar resolvent response of the declared model and source. It is not
a molecular ground-state energy, a ground-state spectral function, a continuum
limit, or a statement about the accuracy of the Hubbard model for a real sample.
Energies are expressed in hopping units; response radii are inverse-energy units.

## Construction mechanism

Local site states are empty, up, down, and doubly occupied. Tensor bonds carry
exact cumulative up/down particle counts. Within each allowed charge block,
QR and singular-value decompositions construct new correlation coordinates.
A finite-state Jordan-Wigner operator encodes the supplied Hamiltonian.

The seed uses truncated complex-symmetric Krylov steps. Local Galerkin equations
then optimize the response tensor in its current bond spaces. When an exact
bound misses the target, the solver forms the physical residual and uses it to
expand the tensor spaces before further local solves. It retains the best
exactly checked attempt.

Operator applications and accepting contractions stream across sites. They do
not construct a complete determinant basis or retain the expanded operator/state
chain. `reference.py` deliberately contains a separate full-sector control;
that control is never imported by the production construction path.

Krylov recurrence residuals, discarded singular values, and local GMRES status
are diagnostics. They are not acceptance certificates. The recorded seed on the
20-site example was poor; local response optimization repaired it. Some local
GMRES calls did not converge to their internal tolerance even in successful
runs. The recomputed exact global residual determines whether the result passes.

## Why the check is valid

Let `A = z I - H`, `y = A^(-1) b`, and let `x` be any proposed complex tensor.
Use the stationary estimator and residual

`g_hat = 2 b^T x - x^T A x`, and `r = b - A x`.

Because the real Hubbard Hamiltonian is symmetric, `A^T = A`. Direct expansion
gives `g - g_hat = r^T A^(-1) r`. Hermiticity of `H` and positive `eta` give
`||A^(-1)|| <= 1/eta`. Therefore

`|g - g_hat| <= ||r||^2 / eta`.

The stationary estimator uses transpose; the residual norm uses conjugate
transpose. Confusing them would invalidate the result. The checker computes
five exact contractions: `x^dagger x`, `x^dagger H x`, `(H x)^dagger H x`,
`x^T x`, and `x^T H x`, plus source amplitudes. No assumption about convergence,
normalization of the proposed tensor, or exact numerical canonicalization is
needed.

Every proposed core is rounded to explicit Gaussian integers divided by `2^30`.
The checker evaluates that rounded program exactly with Python integers or GMP,
and obtains rational centers and radii. Rounding, compression, and incomplete
optimization are all represented in the checked residual. The checker binds
the certificate to the requested model, source, frequencies, and particle
sector. Its software, compiler, Python, and GMP remain part of the trusted
implementation; this is not a formal verification of the executable.

## Measured result

A fresh 20-site open two-leg ladder with `U/t=8`, both hoppings `t=1`, and fixed
`eta/t=1.5` met the requested `0.001` response radius at all three frequencies:

| omega/t | Exactly computed radius, shown approximately | Maximum tensor bond |
|---:|---:|---:|
| 2 | 0.000238726756094269 | 96 |
| -1 | 0.000245167573518541 | 96 |
| 8 | 0.000693634220825039 | 96 |

The first checked bond-48 proposal missed the target with radius
`0.005087971777706870`. Its witness and failure remain in the delivery.
The accepted omega=2 witness has 62,672 nonzero tensor entries. The balanced
sector has 34,134,779,536 determinants; a single complex128 vector over that
sector would require 546,156,472,576 bytes. The construction did not allocate
that vector. The three-query run took 191.98 observed seconds including failed
attempts and independent checks. Parent peak RSS was 286,760,960 bytes; the sum
of recorded parent/child peak RSS values was 573,521,920 bytes. This sum is
conservative bookkeeping, not a simultaneous process-memory trace.

A separate inverse-selection demonstration compared rung couplings 3/4, 1,
and 5/4 at omega=2. Its exact intervals certified 3/4 as the best of those three
choices for maximizing `-Im g`. This is a finite parameter choice, not a newly
invented material or manufacturing recipe. Reusing the original seed is
recorded, and its construction cost is separate from adaptation/checking cost.

A 50-site extension was checked with maximum bonds 48 and 96 and two sweeps per
stage. It remained above target: the final exact radius was approximately
`0.012019191263374638`. This is a failure under that tested budget, not a lower
bound on what a larger representation or different algorithm could achieve.

The separate matrix-free full-sector control met the numerical target on 8 and
12 sites. Its 20-site run was rejected by a 1 GB memory preflight. This is a
representation/resource comparison, not evidence of superiority over mature
DMRG or tensor-network response implementations. Some development runs
occurred concurrently, so timings are not isolated speed benchmarks.

The execution snapshot passed 98 tests, including 11 new package tests. These
cover independent CAR bit-action oracles, rational hopping coefficients,
fermionic signs, lazy operator actions, exact Python/GMP agreement, dense
small-system response containment, request tampering, zero Hamiltonians,
insufficient-rank refusal, and separate-process verification.

## Research relationship

Correction-vector/dynamical DMRG and residual-enriched tensor methods precede
this implementation. Relevant primary references are Eric Jeckelmann,
*Dynamical density-matrix renormalization-group method*, Phys. Rev. B 66,
045114 (2002), arXiv:cond-mat/0203500; and Sergey Dolgov and Dmitry Savostyanov,
*Alternating Minimal Energy Methods for Linear Systems in Higher Dimensions*,
SIAM J. Sci. Comput. 36(5), A2248-A2271 (2014), arXiv:1304.1222.

The advance established here is a new working and exactly accepted capability
inside Spectra: construction before sector expansion, with explicit refusal
when accuracy is insufficient. The tensor parameters are optimized per
instance. This delivery does not establish a newly invented tensor algorithm,
a pretrained transferable theory constructor, recursive learning, or ESI.

`RESULTS.json` contains compact measured results and witness hashes. The full
delivery archive contains executable witnesses, exact receipts, process logs,
negative results, original-source preservation checks, and replay records.
