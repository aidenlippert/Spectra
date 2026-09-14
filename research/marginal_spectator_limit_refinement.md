# Refining the spectator-family limit

The W=0 fixed-family gap is now about three times smaller, using a physical
mixture with **77 sources instead of 85**. Independent exact replay also accepts
slightly stronger lower certificates for W=0 and W=1. The current relaxation's
unrestricted numerical optimum remains unresolved, especially at W=1.

This refinement changes discovery and evidence collection only. The production
ENERGY v13 and FAMILY v9 verifiers, their exact constraints, and their refusal
paths are unchanged.

## Selected exact results

The selected artifacts are
`results/marginal_graded_hubbard8/spectator_hopping/<case>/limit_refined/`.
Exact fractions, provenance and receipt paths are in
`limit_refinement_summary.json` and `limit_refinement_provenance.json` at the
spectator-hopping root. The earlier `combined_summary.json` remains a historical
checkpoint, not the latest W=0 or W=1 bound.

| Target | Million-site open lower/site | Physical upper/site | Current periodic family gap | Previous gap | Mixture sources |
|---|---:|---:|---:|---:|---:|
| W=0 | −0.643062725064887 | −0.610676347051188 | 0.0000679071728421 | 0.000209786662010 | 77 |
| W=1 | −0.661009375948134 | −0.618424482369328 | 0.000511071087230 | 0.000512088956062 | 85 |

These are the same open, half-filled, million-site U=4, t=1, V=1/2 targets.
The lower gains over the previous accepted spectator certificates are about
`1.4659211e-6` and `1.0178688e-6` per site, respectively. The physical uppers were
freshly recomputed, including exact 24-site contraction inside a 160-bit interval
enclosure and the million-site recurrence; they are unchanged.

For W=0, the accepted **periodic** lower is approximately −0.643060225064887,
and the family ceiling is approximately −0.642992317892045. These bracket the
best attainable periodic lower within the matched unrestricted family. For W=1,
the prior family ceiling was retained and independently matched to the new lower.
A family ceiling is not a physical ground-energy upper.

W=+0.1 and W=−0.1 retain their preceding `polished/` checkpoints. No new transfer
claim is made for those cases in this refinement. Both selected new certificates
still strictly exceed the older spin-family ceilings, with fresh exact comparison
at identical projector sources, ratio, ceilings, sparse span and physical target.

## What the numerical experiments established

### Spectral cutting planes

Each lowest eigenvector supplies a supporting plane for the concave minimum
eigenvalue. A bounded LP model used those planes in the existing 85-variable
chart. Both runs made very small improvements, then terminated with HiGHS
`Unknown` status. They were not converged optimization results. Their rounded
lower proposals nevertheless passed the independent full exact replay.

The W=0 run used 68 matrix evaluations and the W=1 run 334; their bounded
iteration histories are retained in `cutting/`. The actual improvement never
exceeded the predicted supporting-plane improvement in those histories. This
numerical consistency check does not establish optimality.

### Small-subspace conic discovery

The next approach compressed the local PSD constraints into low-energy
subspaces and enriched them using full spectral checks. The first W=0 trial
stopped at its explicit dimension-24 cap without an improved selected lower.
The subsequent two runs used at most 12 rounds, 32 active blocks, and compressed
dimension 48. Eight blocks were active in practice. Each solve requested a
15-second solver limit and 60 iterations; these are solver settings, not hard
process wall-clock guarantees.

A temporary coefficient box of radius 0.05 around the seed kept these discovery
models small. Both dimension-48 runs reached their round limit. Their best
rounded candidates became the selected exact lower certificates only after
all 94 physical blocks were rebuilt and checked. The exact replay still spans
4096 states, has maximum local PSD dimension 200, and retains all 151 PSD checks.
The compressed numerical cones did not replace those checks.

The first attempts completed their numerical work but failed when exporting
tuple-keyed bases to JSON. The serialization was fixed, the failed source and
logs were retained, and both bounded runs were repeated successfully. Two
dependent atom-export invocations also failed because the earlier files did
not exist; they were rerun after successful export. These were execution and
serialization errors, not positivity failures.

### Why the tight-looking temporary-box objective was insufficient

A frozen-basis conic solve supplied numerical dual eigenvectors. Eight dominant
vectors were used per target. Their approximately normalized mixtures failed
the full-family zero-moment conditions: the largest numerical residuals were
about `1.85e-5` for W=0 and `2.24e-4` for W=1. The temporary coefficient bounds
matter; a dual for that bounded model need not be a dual for the unrestricted
family. These mixtures were not accepted as family certificates.

The vectors were instead treated as untrusted candidate sources. Adding 200
rounded or deterministically perturbed vectors to the existing sampler and
solving the original moment equations produced a valid W=0 mixture with 77
sources. It uses sampling scale 0.0002, determinant anchors, and **zero pricing
rounds**. Every moment, physical expectation, fidelity bound and source-count
gate passed exact production replay. The corresponding W=1 proposal was weaker
than the previous accepted ceiling and was not selected.

### Wider conic runs and a concrete refutation of incomplete objective values

Removing the temporary box initially caused large full-matrix violations and
reached a total compressed-order cap of 1000. A follow-up added all 4096
projected diagonal nonnegativity inequalities to the numerical model, retained
all 94 initialized blocks, and allowed total compressed order up to 2000 with
maximum cone dimension 48. These are necessary consequences of the same PSD
conditions, not new physical constraints.

Those wider runs reached their 12-round limit without improving the selected
lower. Several final solves returned `user_limit`, with negative full-matrix
slacks. Some reported objectives were even below an already accepted exact
periodic lower. `solver_limit_refutation.json` records exact rational
contradictions for two W=0 values and three W=1 values. Consequently those
displayed numerical values **cannot be global family ceilings**. This is a
refutation of using incomplete solves as optimality evidence, not a failure of
convex duality or a determination of the actual optimum.

The final requested 15-second solver limits were observed at roughly 16.29
seconds and 15.46 seconds in the last rounds. No hard wall-clock or performance
scalability claim is made.

### Exact-basis recovery remained unsuccessful

With 40 pricing rounds and two eigenvectors per priced block, the conic-atom
sampler proposed better numerical ceilings but selected 83- and 84-column
bases that were inconsistent over the rationals. No certificate was accepted
from either run.

Repeating the final LP with equality-row scales 100, 10000 and 1000000 did not
repair the exact basis. The equations themselves were unchanged. Explicitly
adding the omitted fidelity-slack columns also failed: the resulting bases
were dependent, had negative exact weights, or remained inconsistent. The
refusal paths were preserved; no equation was dropped and no negative weight
was clipped.

Both bounded candidate sets are saved in
`scaled_basis/candidate_ledger.json.gz`, including physical integer vectors,
sector coefficients, numerical rows, energies and the exact right-hand side.
They can support a different basis-selection attempt without repeating physical
sampling or spectral pricing. The ledgers and all discovery diagnostics are
nonaccepting evidence.

## Remaining physical obstruction

All 14 spectator moments vanish exactly in both selected witnesses, as do the
previous spin and unconditioned hopping moments. Independent ablation confirms
that deleting only the spectator correction makes the submitted local residual
negative, while restoring it makes the same integer vector's residual
nonnegative.

All four pair-transfer moments remain nonzero. For
\(J(i,j)=d_i^\dagger d_j+d_j^\dagger d_i\),
\(d_i^\dagger=c_{i\uparrow}^\dagger c_{i\downarrow}^\dagger\), the reflected
five-site differences and six-site translated telescopes give:

| Pair | New W=0 witness | Selected W=1 witness |
|---|---:|---:|
| (0,1) | +0.00119448813820 | +0.00101744345132 |
| (0,2) | −0.000204923601442 | +0.000189622342766 |
| (0,3) | +0.000256437762673 | −0.00128950826226 |
| (1,2) | −0.0000269768524881 | −0.000185370967011 |

The W=1 mixture is unchanged, so its moments are unchanged. The exact rank-four
independence argument, full-Fock CAR and symmetry checks, periodic cancellation,
and telescope norm bound of four were rerun. These conditions exclude a
stationary quantum extension matching the averaged six-site density matrix;
they do not exclude every inhomogeneous extension of one marginal. Pair transfer
has not yet been integrated into the energy certificate.

The averaged charge laws still pass the exact stationary classical Markov
extension checks. This does not establish quantum representability or a global
fixed-particle-number extension.

## Validation and continuation

Twenty new final receipts and five intermediate accepting receipts were audited
alongside the previous 40 receipts: **65 receipts and 2200 source-hash entries**
match. No production verifier changed. The prior full regression of 847 tests
and 102 subtests was not rerun; this refinement was verified with fresh exact
physical and family replays plus independent obstruction checks.

No local calculation or GPU remains running at this checkpoint. No GPU was
used in this refinement. Timing observations from concurrent discovery runs
are not controlled benchmarks.

The next numerical work should select a genuinely feasible exact basis from
the saved candidate sets, or improve the bounded conic approach while retaining
full physical verification. The four independent pair-transfer directions
remain an alternative route to stronger compact certificates. The W=1 limit
is still loose. General quantum representability, molecular/long-range/higher-
dimensional transfer, and requested-accuracy scalability remain unproved.
The broader goal remains active.
