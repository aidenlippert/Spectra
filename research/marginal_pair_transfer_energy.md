# Pair-transfer constraints produce strictly stronger energy certificates

Four previously violated pair-transfer consistency conditions now enter the
production energy verifier. Exact replay accepts stronger lower certificates
for W=0 and W=1. They exceed independently replayed ceilings for the entire
preceding fixed spectator family, so the improvement cannot be obtained merely
by optimizing the old correction coefficients more carefully.

Authoritative artifacts are under
`results/marginal_graded_hubbard8/pair_transfer/`. The selected W=0 directory is
`W_zero/`; the selected W=1 directory is `W_plus_1/polished/`. Exact fractions,
current-source receipts and validation details are recorded in `summary.json`.

| Target | New periodic lower/site | Old family ceiling/site | Strict excess above old ceiling |
|---|---:|---:|---:|
| W=0 | −0.642978716055446 | −0.643057945542764 | 0.0000792294873185 |
| W=1 | −0.660850491657565 | −0.660921298140620 | 0.0000708064830552 |

The original targets remain half-filled U=4, t=1, V=1/2 chains. A family ceiling
is an upper limit on the lower certificates obtainable with fixed projector
sources, ratio, ceilings and sparse correction span. It is not a physical
ground-energy upper bound. The new pair-augmented family's own numerical limit
has not yet been certified.

## The added constraints and proof

Define `d_i†=c_i_up† c_i_down†` and `J(i,j)=d_i†d_j+d_j†d_i`. For the pairs
`(0,1), (0,2), (0,3), (1,2)`, define the five-site operator

\[
Y_{ij}=J(i,j)-J(4-j,4-i),
\]

and the six-site difference `T_ij=Y_ij,left−Y_ij,right`. Its sum over translated
windows is exactly zero on the periodic chain. These four real coefficients
are added to the local matrix before the exact positivity check. No physical
pair interaction is added to the target Hamiltonian.

If the corrected local operator plus the nonnegative projector penalties is
at least `ell I` on the full Fock space, summing its translates cancels all
telescopes. The independently proved overlapping-projector ceilings subtract
the same penalty costs as before. Division by five gives the periodic lower
density. Removing the physical closing interactions gives the open-chain
bound with unchanged cost `(2t+|V|+2|W|)/N`.

ENERGY v14 requires a nonzero exact pair-transfer coefficient dictionary and
retains every previous v13 correction and validation gate. Old versions refuse
the new field explicitly. Independent bit-swap tests reconstruct the pair
operator on all 4096 Fock states and verify the fermionic signs, PH and spin-flip
symmetries. Tests also reconstruct the entire physical image from the reflection
blocks and verify translated cancellation. A large pair coefficient is refused
by the PSD test, confirming that these terms enter the accepted matrix.
Fresh probes of the preceding ceiling mixtures also recheck all four nonzero
pair moments, the exact rank-four independence test and the norm bound of four.

The exact energy replays retain all 94 local symmetry blocks, total dimension
4096 and maximum local PSD dimension 200. The number of local blocks and the
overlapping-projector constraints are unchanged. The additional operator data
is four scalar coefficients. This is a finite compact extension; it is not a
scalability theorem for increasing support or requested accuracy.

The existing FAMILY v9 verifier was not enlarged. It now also refuses fixed
pair coefficients, and its matching driver refuses v14 energy certificates.
This prevents an old-family ceiling from being presented as a ceiling for the
larger family. Separate `strict_family_separation.json` receipts check the
unchanged target, chain size, projector data and sparse span, audit the fresh
proof inputs, and compare the lower and old ceiling with rational arithmetic.

## Discovery and exact acceptance

The numerical search used the previous full-spectrum thermal method with four
additional coordinates, for 88 total. Each run had a hard 500 full-spectrum
evaluation budget, two temperatures, and 120 iterations per stage. The search
retained finite-difference gradient/Hessian checks and fresh physical matrix
reconstruction. Numerical success or stationarity was never an acceptance gate.

Initial searches used 334 and 451 evaluations for W=0 and W=1. Both terminated
at iteration limits. Subsequent bounded polishing used 301 and 412 evaluations.
The rounded W=0 proposal worsened and was retained only as a failed improvement
attempt. The W=1 polishing proposal improved the lower certificate and passed
fresh exact replay. No enlarged-family optimum follows from these searches.

The numerical eigenvector changes of basis are nonaccepting auxiliary data.
The standard-library replay reconstructs all physical matrices and checks each
integer congruence inequality exactly, along with the full projector-overlap
proofs. The old family ceilings were freshly replayed under the current code
in the selected directories' `previous_family/` subdirectories.

The selected million-site open-chain energy intervals are:

| Target | Open lower/site | Physical upper/site |
|---|---:|---:|
| W=0 | −0.642981216055446 | −0.610676347051188 |
| W=1 | −0.660854991657565 | −0.618424482369328 |

The physical upper endpoints are unchanged and freshly verified using the
fixed filtered state. Exact 24-site contraction is checked against a 160-bit
outward enclosure before the million-site transfer calculation.

## Frozen transfer to a new coupling combination

A separate experiment changes all three interaction parameters to
`U=5, t=1, V=1/4, W=−1/5`, on the same half-filled chain geometry. The projector
sources, penalties, correction coefficients and trial-state filter are copied
from the initial W=0 pair certificate. Physical profiles are rescaled to the
new interactions and only the scalar local threshold is recomputed. No
coefficient search is performed at the new target.

Independent replay accepts the transferred physical interval:

\[
-0.523070106055446\ \le E_0/N\le -0.488561680298955.
\]

Removing only the pair-transfer terms, leaving all other coefficients fixed
and recomputing the scalar threshold, gives the weaker accepted lower bound
`−0.547307406055446` per site. Keeping the pair terms therefore improves these
fixed-coefficient certificates by exactly `0.0242373` per site. The
`held_out/frozen_transfer_comparison.json` receipt checks both certificates and
the complete frozen-transfer recipe exactly.

This demonstrates transfer of a useful certificate recipe across a new coupling
combination. The comparison does not bound a reoptimized no-pair family. The
geometry, half filling and short interaction range remain unchanged; generic
molecular systems, long-range interactions and higher dimensions were not
tested in this experiment.

## Validation, provenance and remaining work

The full regression passed 898 tests and 102 subtests in 620.84 seconds.
The focused operator/version/refusal regression passed 67 tests. Nine additional
separation-scope and old-family-driver refusal tests passed separately after
full-suite collection. The full log, JUnit XML and `summary.json` record these
completed runs. All numerical and verification processes are terminal.

Seventeen current receipts and 613 source-hash references were audited against
the current files. The provenance manifest preserves the selected proofs,
nonaccepting search trials, source code and validation outputs.

Production source changed in this continuation. The preceding 97 receipts are
historical evidence, not current-source replays. Their 3100 hash references were
audited: 2977 still match current files, and 123 match preserved pre-change
source snapshots in `source_before/`. Fresh selected energy and old-family
receipts establish the current results. Historical files were not overwritten
to conceal source changes.

The pair-transfer energy-integration blocker is resolved for these finite
targets. A compact certified ceiling for the enlarged family remains to be
constructed. Further quantum consistency conditions may still be missing;
there is no general representability oracle or proof of requested-accuracy
scaling. No GPU, paid resource or controlled performance benchmark was used.
The broader goal remains active.
