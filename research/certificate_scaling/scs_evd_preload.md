# SCS `dsyev_` preload experiment

Root replaced the initial unvalidated draft. `scs_evd_preload.c` exports the
LP64 `dsyev_` ABI used by SCS and delegates to `scipy_dsyevd_`. It queries and
allocates both real and integer workspaces, validates arguments and allocation,
and reports separate real-call/query/failure counts. It loads the backend with
`RTLD_LOCAL|RTLD_DEEPBIND` using the absolute `SPECTRA_EVD_LIBRARY` path.

The actual defined symbol (`T scipy_dsyevd_`) is in SciPy's bundled
`libscipy_openblas-c128ec02.so`. The earlier `U` in cython_lapack was only an
undefined reference. No package or installed library is modified.

Two failures were diagnosed before acceptance. The original agent draft had
an undersized integer workspace. Directly linking the corrected preload to
SciPy OpenBLAS then caused a separate symbol collision: gdb traced old SCS
OpenBLAS `idamax_` into an incompatible new-library dispatch routine. Lazy
local loading fixed that collision. Failed logs are retained and do not count
as successful controls.

Build on the existing B host, under `/tmp/spectra-evd-gate`:

```sh
cc -std=c11 -O2 -Wall -Wextra -Werror -fPIC -shared scs_evd_preload.c -ldl -o libspectra_evd.so
```

Set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` before process startup,
`LD_PRELOAD` to the absolute shim path, and `SPECTRA_EVD_LIBRARY` to the
absolute SciPy library path. Never globally preload or replace a library.
The shim is an isolated Linux LP64 research experiment, not a general BLAS
replacement or a supported production solver distribution.

Executed controls (`scs_evd_control.py`): stock and shim SCS each solved 4x4
and 20x20 known PSD problems in 125 iterations. Objectives differed by less
than 3e-14; absolute reference error stayed below 7.05e-8 and minimum cone
eigenvalue above -6.85e-9. The shim reported 504 real calls, two queries and
zero failures. An actual H4 full-cubic solve then passed independent
`python -S` replay with interval width 3.89191924e-7 Ha, 134442 intercepted
real calls and no backend failures. H4 total discovery/export was 6.97 s.

Receipts and failed-debug logs are retained under
`results/certificate_scaling/scs_refinement/preload_gate/`.
No H10 speed or accuracy claim follows from these small controls. Checkpoint
identity now records declared preload/backend content hashes to prevent a
silent restart across a changed experimental backend.
