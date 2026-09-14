# Full number-ideal multiplier ablation

Terminology: `spin_weight2(w)==0` restricts X to Sz-weight zero. It is not a
full SU(2)-invariance test. The SOS Gram blocks are the SU(2)-invariant part;
the two multiplier variants below are called `restricted_Sz_X` and
`full_number_X`.

Provenance: the local H4 and H6 adapter runs below used the repository's
`spin_irrep/*/symmetric_hamiltonian.json` inputs, so they are nearby
symmetrized-H diagnostics rather than original active-space-H endpoints.
Root's clean same-host original-H6 matched run supersedes those preflights;
its artifacts are under
`results/all_angles_20260913/spin_matched/{h6_restricted,h6_full}/`.

`ablation.py` executed a matched H4/H6 basis comparison using the frozen
original-H certificates and the repository's exact canonical maps.  The
physical singlet sector was not substituted into either calculation.

| fixture | full number-conserving degree-2 ideal | SU(2)-invariant subspace | extra full terms |
|---|---:|---:|---:|
| H4 | 235 | 111 | 124 |
| H6 | 1174 | 499 | 675 |

The full ideal is therefore a larger coordinate basis (a strict cone
inclusion is not asserted because number-sector quotient relations may create
redundancy).  This does
not imply a better coefficient-L1 lower bound: coordinatewise L1 is not
SU(2)-invariant, so twirling preserves the zero-residual SOS cone but can
change the numerical objective.  The source multiplier term counts (185/919)
are independent of the raw basis counts because the solved polynomials are
sparse.  Every reported comparison retains the original Hamiltonian and
would be checked by the exact full verifier; no solver status is treated as a
certificate.

Follow-up experiment: run the same `spin_invariant_discovery.run` with only
the `spin_weight2(w)==0` predicate removed from `allowed`, keeping parity,
word groups, row scaling, solver, and time cap identical. Compare exact
replay width and multiplier coefficient L1. A strict improvement would show
that the invariant restriction itself costs objective value; equal values
would suggest the restriction is numerically harmless at that budget. The
larger map and 2.1x/2.35x ideal dimensions should be reported as discovery
cost, even if the objective improves.
## Executed H4 matched solve

Both variants used SCS, `eps=1e-8`, and a 25-second cap. The full-X adapter
changed only the multiplier predicate; Gram decomposition and original-H
verification were identical.

| variant | ideal variables | coefficient rows | wall | exact lower | residual L1 | factor rows |
|---|---:|---:|---:|---:|---:|---:|
| invariant X | 111 | 777 | 8.88 s | -3.6670007606 | 8.1663e-7 | 705 |
| full X | 235 | 1,809 | 7.94 s | -3.6670008024 | 8.5424e-7 | 679 |

At this matched budget full X is marginally weaker after exact export (about
4e-8 Ha) and has a slightly larger residual. This is one H4 sample, but it
shows unrestricted number-ideal multipliers do not automatically cure the
invariant solver while increasing affine-map rows by 2.33x.

## H6 execution

The configured SSH aliases `GPU_math` and `lambda` were unresolved from this
workspace, so HOST B could not be reached. I ran the full-number-X adapter
locally with the requested 120-second SCS cap. It built 1,174 multiplier
variables and 20,021 coefficient rows, reached `optimal_inaccurate` at the
cap, and exact original-H replay gave lower `-6.3332214618` Ha with residual
L1 `1.4385e-4` Ha. Wall time was 139.50 s including export. This is directly
comparable to the existing restricted-Sz H6 run only as a diagnostic because
that archived run used a different cap; a matched restricted-Sz rerun is the
next needed control once a host is available.

The clean matched original-H6 comparison subsequently completed by root gives
restricted lower `-6.333242734124332` Ha and full-number lower
`-6.3332400402471665` Ha, an improvement of only `2.693877e-6` Ha for the
larger multiplier basis. Those results supersede the local symmetrized-H
preflight for endpoint claims.
