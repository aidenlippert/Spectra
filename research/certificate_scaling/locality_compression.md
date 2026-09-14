# Locality-compression certificate probe

This is a bounded falsification experiment for the proposed structural idea:
small spatial clusters might produce useful global lower certificates.  The
Hamiltonian is the open, spinful one-dimensional Hubbard chain

\[
 H=-t\sum_{i,\sigma}(c^\dagger_{i\sigma}c_{i+1,\sigma}+h.c.)
   +U\sum_i n_{i\uparrow}n_{i\downarrow},
\]

at half filling.  This is a deliberately clean fermionic family; it is not a
molecular basis-set calculation.  The script solves each fixed-particle
sector by a numerical sparse eigensolver (so cluster values are numerical
lower-bound candidates, not independently rounded-proof certificates).  For a partition into clusters of
maximum size `k`, it minimizes the sum of independent cluster ground energies
over all particle allocations.  If there are `b` cluster boundaries, the
omitted hopping operator is bounded by `2 |t| b`: one unit operator norm per
spin and boundary.  Therefore

`lower = sum(cluster energies) - 2 |t| b`

is the exact mathematical form once cluster minima are certified; here only
the coupling penalty is exact.  The same decoupled cluster state is a
variational product upper bound, so the reported omitted-penalty width is
`2 |t| b` subject to the numerical cluster qualification. Every state
count, sector minimization, penalty, and wall time is recorded in
`summary.json`; no full-sector calculation is hidden inside the certificate.

The script also compares overlapping windows. Windows start every `k-1`
sites, and each onsite and bond term is assigned to exactly one owner window;
the resulting local Hamiltonians sum to the global Hamiltonian, so the sum of
their independent minima is a lower-bound candidate with no explicit cut
penalty. Since overlapping windows do not share one particle sector, local
particle numbers are minimized independently. This is a stress test rather
than a claim of optimality.

## Result

The sweep used `L=4,6,8`, `U/t=0,2,4`, and cluster budgets `k=2,3,4`.
The gaps `E_exact - lower` for `U/t=2` were:

| L | k=2 | k=3 | k=4 |
|---:|---:|---:|---:|
| 4 | 1.5962 | 0.9441 | 0 (full cluster) |
| 6 | 3.1619 | 1.0939 | 1.5657 |
| 8 | 4.7186 | 2.6506 | 1.5263 |

For fixed cluster size, the gap grows roughly with the number of boundaries.
The omitted-coupling norm penalty dominates the hoped-for locality gain.  Even
when `k=4` is fixed, the certificate is not approaching a size-independent
absolute error over this range.  The `k=3` irregularity is a finite-size
partition effect, not evidence of a theorem.  At `k=L` the lower bound is the
exact diagonalization result by construction.

This falsifies the naive claim “short-range clusters plus a triangle penalty
scale automatically.”  It does **not** rule out a stronger structural result:
one would need a boundary correction that exploits state-dependent density
matrices, screening, cancellations, or a small interface certificate rather
than paying the operator norm of every omitted hopping.

Run the smoke sweep with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python research/certificate_scaling/locality_compression.py \
  --sites 4 6 8 --clusters 2 3 4 --couplings 0 2 4
```

It produced 27 records in about 1.4 seconds on the local CPU.  A larger run
with `L=16,32,64`, `k=2,4`, and `U/t=0,2,4` is supported without global
diagonalization using `--global-max-sites 8`; it reports product upper bounds
and cluster lower candidates only.  The script accepts `--backend cupy` for
dense cluster sectors, failing clearly if CuPy is unavailable.  This supports
hundreds of small cluster parameter cases, not a claim that GPU execution
makes the global problem scalable.

## Scope and failure modes

The chain is one-dimensional and finite, and the exact diagonalization is only
an auditing reference.  This does not establish an asymptotic lower bound for
molecules, nor does it test orbital locality in a Gaussian basis.  It does
establish an explicit baseline against which a proposed compressed interface
certificate must improve: same Hamiltonian, same cluster budget, and a penalty
whose omitted operator norm is fully accounted for.
