# H8 fixed-number completion: measured outcome

The best new complete interval is **2.849197 mHa**. It does not meet the 1.6 mHa target. Its change relative to the retained paired interval is +0.025759 mHa of narrowing.

The supplied reduction is now implemented and checked against the actual strong
H8 proof. The new constructor permits independent cubic Gram blocks and collective
fixed-number cancellation. This is a completed test of the proposed algebraic
change; it is not a general many-body solution or a proof that every compact
representation has been exhausted.

## The actual strong proof confirms the missing freedom

The inherited proof contains 2,848 quartic number-multiplier coefficients. Exact
CAR expansion verifies **W+L2(X4)=-R6**, with the project's actual sign convention.
Projection verifies the corresponding V and contraction-free Z relationships
exactly. Coefficient Frobenius norms are 0.4608828142 for W, 0.1318644498 for V,
and 0.0000141065 for Z. The last quantity is nonzero: it remains covered by the
old exact residual witness, never silently discarded. These are coefficient norms,
not molecular excitation gaps or energy improvements.

The independently replayed inherited full proof remains at
**1.176860 mHa** using the frozen upper. It was read for diagnosis
only; no full-proof factors or eigendirections entered the new compact spans.

## Complete H8 comparisons

| Construction | Full width, mHa | At most 1.6 mHa? |
|---|---:|---|
| Preserved strong full-cubic proof | 1.176860 | Yes |
| Preserved compact paired proof | 2.874956 | No |
| Restored freedom, fixed 1.4 mHa proposal | 3.558913 | No |
| Restored freedom + prescribed linear completion | 2.849197 | No |
| Restored freedom, original compact spans | 2.896071 | No |

Every new row has a separate standard-library exact replay using the sealed
spin/CAR verifier. It includes the S=0 proof, the separate MS=1 proof covering
all nonsinglets, and the original-H defect 59/250000000000 Ha charged once.
All rows use exactly the same frozen rational upper endpoint. The same valid
upper is reused here; its prior construction and verification are dependencies.

The best new candidate has a raw exported interval of 2.843009 mHa;
the complete residual and original-H spin adjustment adds 0.006188
mHa. Its remaining distance from the target therefore cannot be attributed
mainly to replay precision. The final improvement belongs to the combined change
with the prescribed linear directions. A matched paired run with those appended
directions was not performed, so the energy gain cannot be assigned uniquely to
unpairing rather than the additional directions or further optimization.

The strict unchanged-span model has 91,752 independent Gram entries, plus
18,128 for the nonsinglet proof. The paired comparison has 54,932 plus 18,128.
The old paired construction is a feasible special case. Its numerical embedding
matches every enlarged coefficient equation to 4.29e-10. All retained cross terms
remain available. No higher MPS moments were filled in or treated as known.

After the unchanged-span test, an additional fixed-number issue was identified:
linear words had been removed using identities available in the full cubic space,
but their replacements were 94.83–97.59% outside the compact spans in relative
Euclidean distance. An exact completion appended 16 prescribed number-dressed
linear directions and their adjoints. This increases the singlet optimization to
96,904 entries, or 115,032 including the nonsinglet proof. It retains every old
column. The additions follow the number identity; no larger-rank ranking search
or full-proof teacher chooses them.

## What the numerical searches establish

The direct 1.4 mHa target did not survive its residual allowance: its complete
replay is 3.558913 mHa. A number-projection repair retained all contraction-free
remainder and used a separately verified residual wedge witness. Its raw target
was never presented as a certified energy interval.

Early optimization runs had poorly conditioned linear solves. A measured
preconditioner reduced the fixed-target solves from repeatedly hitting 120
iterations to usually converging in one. Nevertheless, its coefficient residual
leveled off. Additional bound-maximization and the prescribed linear completion
were tested. The table contains only accepted full intervals; other saved
iterates remain untrusted proposals, even when a numerical process exited normally.

There is **no exact dual obstruction**. A diagnostic dual had negative moment
eigenvalues and an unbounded-for-this-purpose functional coordinate around 29.77;
it fails the accepting conditions for a residual-compatible family ceiling. The
finite runs therefore do not prove that these spans, much less all compact
representations, can never reach the target.

The projected lift coefficient provides an additional diagnostic of whether the
search actually uses the restored freedom:

| Construction | Coefficient norm of V | Coefficient norm of Z |
|---|---:|---:|
| Strong inherited proof | 0.1318644498 | 0.0000141065 |
| Restored freedom + prescribed linear completion | 0.00126808547 | 5.34881073e-07 |
| Restored freedom, original compact spans | 2.22758527e-06 | 9.35801286e-07 |

These are exact rational squared norms, displayed after taking a floating square
root. They are not energy errors or an objective to match: a different successful
certificate need not resemble the strong proof. Every nonzero Z is paid by the
complete replay. The norm comparison is only a diagnostic of the returned factors.

In the original compact spans, the returned V norm is about 0.00000223. After
adding the prescribed number-derived linear directions, it rises to about
0.00126809. Thus the additional directions do let the returned construction use
more of this freedom. The complete energy interval, rather than that increase,
determines whether the change is useful.

## Costs and dependencies

- 28 completed instrumented processes total **3165.034 process-wall seconds**,
  with 3104.419 recorded CPU seconds. This includes failed and unsuccessful
  attempts, preparations, tests, repairs, and exact replays. It excludes editing,
  reading, interactive diagnostic probes, final inventory, and prior discovery;
  it is not a fresh-problem end-to-end timing claim.
- Peak recorded process memory is **1.483 GB**.
  Each dense normal-matrix cache occupies 582,496,712 bytes. It is a numerical
  construction cost on coefficient equations, not a determinant matrix or a
  small accepting certificate. Two different span caches are retained.
- Fresh full-degree operator preparation took 178.657 process seconds. It
  visits 1,206,592 cubic dictionary word pairs and stores 2,009,792 sparse-map
  nonzeros. The exact spin reduction retains 8,533 independent equations from
  28,461 full coefficient rows. There is no N-particle determinant enumeration.
- The compact-span seed inherits earlier H8 searches totaling 380.135 internal
  seconds, plus direct quartic-guide preparation and the nonsinglet construction.
  The upper MPS inherits its 414.36-second recorded discovery path. Historical
  preparation and earlier failed attempts are not made free by using frozen files.
- The old full-cubic search's 721.563-second historical discovery is a dependency
  of the inherited strong baseline, not of the new constructor. Its diagnostic
  replay in this pass is separately charged. These timings are not a matched
  hardware speed comparison.

The first quotient repair took 177.721 seconds. Proposal-side spin projection
was then accelerated using exact canonical spin flips and tested against the
original projector. The accepting verifier was left unchanged.

## Validation and preserved milestones

The focused suite passes **15 tests**, including all 325 pair-matrix units at
m=5,6, sparse m=12,16 cases, the nonzero zero-contraction counterexample, the
signed-factor norm test, independent six-mode bit action, a nonzero completion
through the unchanged accepting checker, all molecular adjoint-pair maps, all
662 six-mode matrix units for spin projection, and the number-dressed linear
identities. The 64 enumerated bit inputs belong only to an isolated algebra test.

All **11,077 inherited files** and both supplied request files match
their original hashes. The original H6 result and fresh-geometry 1.034056 mHa
milestone remain preserved; the latter is 1.034064 mHa when its actual MPS upper
is independently rechecked with outward rounding. This pass uses the same frozen
H8 upper as its paired comparison.

The mathematical construction and its exact scope are in [MATHEMATICS.md](/Users/aidenlippert/Documents/Spectra/research/sector_quotient_20260914/MATHEMATICS.md).
Machine-readable results are in [audit.json](/Users/aidenlippert/Documents/Spectra/results/sector_quotient_20260914/audit.json), with all process
receipts retained in the runs directory. The [best new certificate](/Users/aidenlippert/Documents/Spectra/results/sector_quotient_20260914/candidates/linear_closure_repaired/certificate.json)
has a separate [complete replay](/Users/aidenlippert/Documents/Spectra/results/sector_quotient_20260914/candidates/linear_closure_repaired/interval.json).

The useful result of this pass is a tested constructor that restores the missing
fixed-number freedom, together with a measured molecular outcome. Further claims
of compact H8 accuracy or a limitation theorem must follow accepting certificates.
