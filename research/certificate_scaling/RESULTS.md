# Certificate compression: first parallel experiment campaign

12 September 2026. Eight subagents and the coordinating agent investigated
operator selection, factor compression, locality, residual domination, constraint
repair, structural assumptions, literature, and verification. This is a completed
bounded experiment campaign, not a solution to general certificate discovery.

## Verified result

Rounding the existing sparse SOS factors to a common denominator of 100,000,
removing zero rows and unused word columns, produced smaller complete molecular
certificates. Each retained factor remains an exact rational square. The existing
CAR checker charges the entire new coefficient residual. A separate checker
recomputes the variational upper endpoint from integer amplitudes and the exact
Hamiltonian. Compression never uses those amplitudes.

| Finite-basis model | Original interval certificate | Compressed interval certificate | Size reduction | Certified interval width |
|---|---:|---:|---:|---:|
| H4 square, side 2 Å, STO-3G | 175,912 bytes | 119,096 bytes | 32.30% | 0.00127285521019 Ha |
| H4 rectangle, 1 by 1.5 Å, STO-3G | 129,032 bytes | 112,873 bytes | 12.52% | 0.000731085390761 Ha |

Both pass the explicit 0.0015 Ha target. Denominator 1,000,000 keeps tighter
intervals: 0.000151990404112 Ha and 0.0000771861367607 Ha, with smaller size
reductions. These results concern rational finite-basis electronic Hamiltonians.
Continuum/basis error and the underlying numerical integral algorithm are outside
the interval. They do not demonstrate new molecular discovery or a new theorem.

The 100,000 denominator was selected from the square experiment before its
unchanged application to the rectangle. The broader precision sweep is also
retained for both. This is limited transfer of a post-processing rule, not a test
of efficient discovery on a larger active space. Both cases have eight spin
orbitals and four electrons. Existing certificate discovery costs remain part of
the method; upper validation explicitly enumerates 70 fixed-N states.

Authoritative exact intervals and complete certificate files are in
`results/certificate_scaling/intervals/summary.json`. Reproduce their replay:

```sh
python -S research/certificate_scaling/certified_intervals.py
python -S -m experiments.marginal_transfer_verify results/certificate_scaling/intervals/square_d100000_interval_certificate.json
python -S -m experiments.marginal_transfer_verify results/certificate_scaling/intervals/rectangle_d100000_interval_certificate.json
```

## What the other attacks established

**Spectral rank truncation is an inadequate default on this fixture.** Eight
rank caps were tested on each of two representations of the square H4 model,
with a CPU comparison for one representation. Every emitted full certificate
was checked exactly. On the GPU-produced unsplit square certificate, rank 32
gave a lower endpoint near -3.617703 Ha, versus -3.330389 Ha for the full factor
set. Smaller matrix rank alone does not establish a useful size/accuracy tradeoff.
Degenerate eigenspaces can also rotate differently across backends; rounded
factors need not be identical. Exact replay is authoritative for each artifact.

**Constraint repair recovers some lost accuracy at fixed factors.** A sparse
linear program re-optimizes the number-conserving multiplier and scalar shift.
Its proposal is rounded and replayed exactly. Lambda results:

| Rank cap | Original compressed lower | Repaired lower |
|---|---:|---:|
| 8 | -10.15146578 | -8.73657976 |
| 32 | -3.61770304 | -3.49463512 |
| 64 | -3.35605340 | -3.34898047 |

These remain worse than the full certificate. Local CPU experiments also
preserve a rounding-induced regression on the full-factor case; the method
must retain the original certificate as a fallback. All quantities above are
lower endpoints, not ground-state predictions. The constraint repair is a
useful representation test, not evidence of efficient global discovery.

**Structured residual control is valid but insufficient here.** The standalone
research checker reconstructs the complete exact residual
`R = H - bI - SOS - (Nhat-N)X`. For `R=cI+dGamma(A)+T`, it centers A by its trace
and uses a rational row-norm bound on the one-body part and coefficient-l1 on
every higher-body word. This improves some unrepaired compressed certificates;
after ideal repair the additional improvement in this batch is negligible.
The production checker is unchanged. See `full_residual_replay.py`.

**Naive locality has an extensive error penalty.** The numerical Lambda sweep
tested 60 settings: lengths 8, 16, 32, 64; cluster sizes 2, 3, 4; and U/t values
0, 1, 2, 4, 8. Each has disjoint and overlapping-window records. The disjoint
method charges `2|t|` for each cut. Its error allowance grows with system size
at fixed cluster size. Overlapping ownership removes this explicit cut penalty,
but independent window minimization can still give loose bounds. These GPU
records are numerical candidates, not exact certificates. Separately, 16 local
integer-coupling cases received rational LDL checks and exact disjoint dynamic
programming composition, through length 64. Their accompanying overlapping
records remain numerical. No global sector is enumerated for the large chains.

**Adaptive selection has not beaten the simple baseline.** On a four-mode
coefficient-space test, all feasible full/fixed/random/adaptive proposals passed
exact export. The fixed three-block dictionary achieved width about 3.81e-7,
versus roughly 1.68e-6 to 2.26e-6 for the adaptive cases. The adaptive search
scans all remaining blocks and pays for many SDP solves. The six-mode Lambda
adaptive run also exported a valid lower certificate. Neither result establishes
a search-cost advantage. The source `hopping_model` is a structured interacting
test Hamiltonian; it must not be confused with an ab-initio Coulomb molecule.

## Compute, validation, and corrections

One Lambda A100-SXM4-40GB was used at the displayed $1.99/hour. Ten bounded jobs
completed in a 42.85-second batch, including one CPU discovery job overlapping
the GPU work. The single end-to-end H4 runs took 5.57 seconds with CuPy and
4.67 seconds with NumPy. This is not a repeated controlled performance benchmark;
it provides no evidence of an end-to-end GPU speedup at this size. Exact replay
and coefficient work remain on CPU.

All 69 downloaded result files matched remote hashes. All 37 returned lower
certificates passed fresh `python -S` replay locally. Six complete molecular
intervals were independently replayed, including the two compressed examples
above. Source and result archives are preserved under
`results/lambda_runs/certificate_scaling/`; `independent_replay.json` and
`lifecycle.json` record validation and resource state. The provider dashboard
was checked after termination and showed **No running instances**. The observed
running interval cost approximately $0.29; this is an estimate, not an invoice,
and excludes any unobserved interval before the first running-state observation.

Several early exploratory artifacts were superseded: a missing dagger in a
per-block CAR map, incorrect LP solution slicing, and a rational-string
reporting error. The complete returned certificates passed the unchanged exact
checker after these corrections. Initial per-block gain claims and false
"export rejection" claims are not evidence. The theorem draft also initially
used an invalid arbitrary two-sided particle-number ideal; it now requires
explicit sector vanishing and includes a counterexample to that shortcut.

The raw remote folder named `heldout_gpu` contains the symmetry-split,
coefficient-rounded representation of the **same square geometry**. It is a
representation comparison, not held-out molecular transfer. Only the subsequent
rectangle experiment supports the limited transfer statement above.

No existing production verifier or solver was changed in this campaign. The
older interrupted v20/v16 chain work and its outstanding full-suite check remain
separate and unfinished; this campaign does not silently complete them.

## Next decisive target

The current best lead is **direct sparse-factor discovery with explicit residual
control**, combining factor support/precision selection with number-constraint
repair. Post-hoc compression proves that some bytes can be removed; it does not
avoid constructing the original large certificate.

The next comparison should enforce a fixed total interval width and count all
candidate generation, SDP/LP work, exact replay, and witness costs on increasing
active spaces. The structural result we seek must bound the cost of finding the
retained factors and the possible improvement from omitted blocks without
enumerating the full omitted dictionary. Neither this campaign nor the
elementary local-norm theorem supplies that result. Broad claims that the
published methods cannot scale are also not established by these experiments.
