# Eight-site Hubbard accuracy target reached

The fully coupled open `2 x 4` Hubbard ladder at `U/t=8` now has a cold,
exactly replayed energy interval of width **4.774093164360818e-7 t**,
below the frozen **0.001 t** target. This covers the ground energy at
fixed total particle number eight, with all spin sectors included by
spin-rotation symmetry.

The result closes this finite accuracy target. The constructor retains all
4,900 amplitude entries, so the requested scalable compression breakthrough
remains open. This fixture is a Hubbard lattice, not the molecular H8 case;
its energies are in hopping units `t`, not hartrees.

## Mathematical change

An exact particle-hole transformation turns the state into a `70 x 70`
amplitude matrix. Spin-reflection positivity proves that a ground
eigenmatrix `W` can be chosen positive semidefinite. For a rational trial
matrix `C`, the two exact conditions

\[
C\succ0,\qquad\mathcal L(C)-\ell C\succeq0
\]

then imply `e_0>=ell` through

\[
0\le\operatorname{Tr}[W(\mathcal L(C)-\ell C)]
=(e_0-\ell)\operatorname{Tr}(WC).
\]

The Rayleigh quotient of the same `C` supplies the upper endpoint. The
proof retains the full amplitude information instead of replacing local
excitations by the coarse bound that limited the preceding pass.

[DERIVATION.md](DERIVATION.md) gives the complete transformation, sign
convention, elementary positivity proof, rational checker and scope.
The reflection-positivity ingredient is established Hubbard mathematics;
see [Lieb (1989)](https://doi.org/10.1103/PhysRevLett.62.1201) and the
[prior-art note](PRIOR_ART.md). No mathematical novelty claim is made.

## Exact result

| Quantity | Accepted value |
| --- | --- |
| Lower endpoint `L/t` | `-30259232831/10000000000` |
| Upper endpoint `U/t` | `-151296140284988781488/50000000000150236249` |
| Interval width divided by `t` | `238704658218758147090919/500000000001502362490000000000` |
| Decimal width | `0.0000004774093164360818` |
| Frozen target | `0.001` |
| Largest accepting PSD matrix | `70 x 70` |
| Amplitude information retained | all `4,900` entries |
| Candidate file size | `43,074` bytes; shared model and checker additional |

The previous local-response interval was `1.0814371852093165 t`. Its
limitation is not contradicted: the new proof belongs to a different
family and retains the full state information. Comparing those two widths
does not establish a fair scaling or runtime advantage.

## Discovery and verification evidence

The cold solve started from the identity matrix. It loaded no saved ground
state or earlier energy endpoint, required 98 amplitude-map applications,
and passed the exact gate on its first rounded candidate. The proposed
floating-point eigenvalue is not trusted by the accepting checker.

The checker reconstructs the model, verifies exact positive definiteness
of the integer `C`, verifies exact PSD of `q L(C)-p C`, and computes the
upper by exact integer contractions. A fresh process with site packages
disabled replayed the proof using only Python's standard library. It
rejected wrong interaction strength, wrong charge populations, an
asymmetric amplitude matrix, a singular amplitude matrix and an
overclaimed lower endpoint.

Five focused tests passed. The independent model-binding test compared
every signed transition to direct physical Fock action on both the
four-site and eight-site models. It enumerated 36 and 4,900 configurations
respectively; that validation cost is included below. The production
constructor and checker do not assemble the `4900 x 4900` Hamiltonian,
but both still use the complete `70 x 70` amplitude matrix. Replay is
independent of the numerical proposer, not an independently developed
second exact-PSD implementation or an external audit.

## Measured cost

All three bounded processes ran on the same local Apple Silicon host,
macOS 26.4.1, with numerical thread counts restricted to one. No external
compute was used.

| Metered process | Wall seconds | Peak child RSS bytes | Outcome |
| --- | ---: | ---: | --- |
| Cold construction, including exact acceptance | 1.229685 | 63,455,232 | Passed |
| Five focused tests, including independent Fock comparison | 2.005998 | 65,191,936 | Passed |
| Fresh exact replay, including five corruption checks | 1.265478 | 15,925,248 | Passed |

The sum of these subprocess times is **4.501161 seconds**. This is neither
the elapsed research-session time nor a new aggregate single-run benchmark.
The constructor's internal timer was 0.789131 seconds. Reasoning, source
inspection, document preparation and earlier failed campaigns are not
contained in these subprocess numbers. Earlier costs remain in their
original receipts. There were no failed bounded runs or rejected candidate
attempts in this positive-cone experiment.

These timings describe a small Hubbard fixture with integer model inputs.
They are not molecular integral-to-certificate timings and do not establish
an advantage over an optimized exact-diagonalization reference.

## Reproduction

From the repository root, choose a new output directory for a fresh
discovery run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv-correlated/bin/python -B -m research.positive_cone_20260916.discover \
  --out /tmp/spectra-positive-cone-new-run
```

The accepted certificate can be replayed without NumPy, SciPy or a solver:

```sh
python3 -B -S -m research.positive_cone_20260916.replay \
  --certificate results/positive_cone_20260916/cold/certificate.json \
  --out /tmp/spectra-positive-cone-new-replay
```

Existing output directories are refused so prior receipts cannot be
silently replaced. The exact interpreter paths and metered commands are
recorded in `results/positive_cone_20260916/runs/*.json`.

The actual proof is
[`cold/certificate.json`](../../results/positive_cone_20260916/cold/certificate.json);
the fresh acceptance record is
[`fresh_replay/receipt.json`](../../results/positive_cone_20260916/fresh_replay/receipt.json).

## Remaining open target

The finite accuracy result is complete. For `s` sites this representation
still needs `binom(s,s/2)^2` amplitudes. A scalable result requires a
directly constructible representation and a rigorous residual-positivity
check that avoid that expansion. The present proof does not establish
those properties or generalize them to arbitrary molecular Hamiltonians.
