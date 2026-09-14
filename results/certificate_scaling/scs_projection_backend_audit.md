# Remote SCS projection backend audit

Both warm hosts use the same installed backend in
`/home/ubuntu/spectra-venv`:

```text
SCS 3.2.8
scs/_scs_indirect.cpython-310-x86_64-linux-gnu.so
libopenblas-r0-11edc3fa.3.15.so (OpenBLAS 0.3.15, DYNAMIC_ARCH, NO_CH, NO_AFFINITY)
```

`ldd` shows the SCS extension linked to the bundled OpenBLAS and gfortran
libraries, rather than the host's system BLAS. `nm -D` shows an unresolved
`dsyev_` reference in the SCS extension. The bundled OpenBLAS exports
`dsyev_`, `dsyevd_`, and `dsyevr_`. This establishes that the installed SCS
extension calls a LAPACK symmetric eigensolver through `dsyev`; it does not
establish that `dsyev` dominates every measured cone-projection cycle.

The OpenBLAS binary identifies itself as version 0.3.15 and includes dynamic
architecture dispatch. Its strings also contain fallback warnings for AVX,
AVX2, and AVX512 kernels. Those strings are compiled diagnostics, not proof
that a fallback was selected on the host. No library loading or symbol timing
was performed while the H10 jobs were running.

The observed H10 stage-0 receipt reports 175 SCS iterations, 666,722.7 ms
solve time, 635,576.4 ms cone time (95.3%), 23,130.2 ms linear-system time,
and 1,480.7 ms acceleration time. The isolated Torch benchmark found an A10
advantage for batch 8×725 eigensolves, but that is a separate kernel and does
not prove SCS can use it.

The smallest useful CPU experiment is an isolated projection microbenchmark
using the exact eight 725×725 and eight 225×225 PSD block shapes and the same
float64 matrices from a saved stage. Compare the current SCS projection path
against direct LAPACK `dsyev`, `dsyevd`, and `dsyevr` calls with one and several
OpenBLAS threads, recording residuals and wall time. Run it outside the live
campaign. If direct calls reproduce the cone time, a narrowly scoped SCS
projection backend change can be evaluated; if they do not, inspect cone
assembly/scaling before changing eigensolver libraries. No conclusion about
unoptimized BLAS or a slow LAPACK algorithm is justified from the current
linkage evidence alone.
