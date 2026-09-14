# Further enrichment improves H6, but the compactness tradeoff weakens

**The certified interval narrowed from 7.414219 to 3.942479 mHa**, a 46.83%
reduction. The continuation used **214 coupled combinations**, up from 64.
The same 214 combinations used separately give **6.619637 mHa**, so coupling
improves the found bound by another **2.677158 mHa**.

The **1.6 mHa target remains unmet**. At 214 combinations, the selected set
is 43.50% of the complete 492-generator frame. This is progress in the bound,
but it does not establish a small sufficient representation of the difficult
correlations. The attempted collective dual repair did not produce an
accepted full-family witness, so the enlarged family's ultimate limit also
remains unknown.

The ten spatial Hamiltonian patterns and their **0.108331264 mHa** certified
collective remainder are unchanged. The growing objects here are the joint
correlation combinations attached to those patterns, not additional
Hamiltonian terms or many-electron wavefunction amplitudes.

## What ran

The six-minute continuation imported the existing CAR maps, pricing,
conditioning, exporter and exact checker without modifying them. The saved
64-direction certificate was independently replayed before its lower was
used to score any gain. Its directions and frame references were required
to agree exactly with the certificate. Its numerical dual supplied proposals
only.

At each round, candidate bundles added one or two independent negative
directions per symmetry block. Both were optimized and exactly replayed when
the remaining budget permitted. Continuation used accepted energy gain per
trial wall time plus pricing cost. All unselected work remained charged.

| Selected combinations | Certified interval, mHa |
|---|---:|
| 64, inherited | 7.414219 |
| 80 | 6.724615 |
| 96 | 6.193385 |
| 112 | 5.679481 |
| 128 | 5.240350 |
| 144 | 4.982903 |
| 160 | 4.586025 |
| 176 | 4.453751 |
| 192 | 4.186035 |
| 206 | 4.087226 |
| **214, coupled** | **3.942479** |
| Same 214, separate contributions | 6.619637 |

Nineteen new candidate trials completed and passed exact acceptance. Nine
additional candidate rows beyond the selected path above are retained in
the ledger. All nine completed two-bundle comparisons selected the larger
bundle; the final round ran only the smaller bundle because of the time
reserve. There is no demonstrated advantage over a fixed larger-bundle
policy. Some blocks eventually supplied fewer than two new independent
directions, which explains the step from 192 to 206 rather than 208.

The search stopped at **332.067 seconds** because only 28 seconds remained
and the declared reserve was 60 seconds. The twelve-round and 256-direction
caps were not reached. Every accepted lower was retained, including a weaker
168-direction candidate that did not replace the 160-direction result.
Finite solver and export accuracy can make a found certificate weaker even
when its mathematical cone contains the earlier cone.

The best molecular lower is -6.337001104797967... Ha. The rational reference
upper is -6.333058626233001... Ha. The best numerical proposal was
`optimal_inaccurate`; its exact residual penalty, **0.004052841 mHa**, is
included in the reported lower. Numerical warnings and statuses remain in
the logs. The rational reconstructed proof, not solver status, accepts it.

## What the extra accuracy cost

| Quantity | Inherited 64 | Continued 214 |
|---|---:|---:|
| Added block dimensions | eight of size 8 | two of size 26, six of size 27 |
| Gram entries including quadratic baseline | 6,236 | 11,450 |
| Numerical map nonzeros | 251,268 | 2,537,322 |
| Compact lower certificate | 94,565 B | 237,366 B |
| Expanded SOS certificate | 645,205 B | 1,453,772 B |
| Nonzero direction coefficients | 3,936 | 13,224 |
| Nonzero compact factor coefficients | 4,736 | 9,402 |

The full frame has 43,920 Gram entries, but its complete solve did not finish
in either earlier control. That does not establish a speedup or an optimum
for the selected representation. The continued proof uses fewer Gram entries,
while its numerical map has become about as large as the earlier full-frame
map without conditioning. The quadratic baseline is still present.

The last selected trial took **48.413 seconds**: 1.227 for construction,
39.622 for solving, 0.039 for export and 7.525 for exact acceptance. The
64-direction trial in the previous pass took about 8.196 seconds. Increasing
the selected set is no longer cheap in this implementation.

The new continuation cost **332.067 seconds**. Including its inherited
descriptor discovery gives **554.578 seconds of cumulative causal discovery**.
The new separate-contribution control cost **24.713 seconds** and inherits
that entire direction-discovery cost. Both together took 356.782 seconds.
Earlier full-family timeouts, diagnostics and validations remain in their
original ledger; none are erased or treated as free successes.

Fresh interval replay took **141.928 seconds**. The separate dual diagnostic
took **11.802 seconds** and failed to produce a witness. These two processes
overlapped; their wall costs are not a sum of elapsed calendar time. All
measurements are single runs, with BLAS thread environment variables set to
one and Clarabel retaining its automatic thread setting.

## The collective repair was uninformative

A separately declared diagnostic tried to repair the final small solve's
moment assignment so that all 492 generators, including every unselected
direction, would have positive moment matrices. This required no full-frame
primal optimization.

After exact affine projection, generalized eigenvalues against the uniform
fixed-N trace proposed a minimum trace mixture of about **1.088947%**.
Every scheduled safety-margin proposal exceeded the declared 1% mixture cap,
so **none entered the complete exact feasibility checker**. No exact dual
ceiling or new error floor was established. The refusal is a diagnostic
outcome, not evidence that the full operator family is insufficient.

Simply raising that cap would not provide the desired accuracy obstruction
through this repair. The first proposed mixture's algebraic objective is
-6.321448901515493... Ha, **11.609725 mHa above the known upper**. That is an
unaccepted objective, not a certified ceiling; even if feasible, it would
not exclude an accurate lower. The uniform-trace repair loses too much
energy information here.

## Best next experiment

The diagnostic suggests testing **trace-normalized direction selection from
the same 64-direction seed under a matched total budget**. The current rule
uses ordinary eigenvectors of moment matrices, normalizing directions by
their frame coefficients. A uniform-trace metric instead normalizes by the
operator's mean squared size in the fixed-N sector. This can expose
violations that are small in coefficient coordinates but large relative to
the operator's size.

This is a testable selection hypothesis, not an established improvement or
a claim that a different metric solves the many-body problem. Its value
must be judged by actual certified energy gain, proof size and complete
construction/solve/replay cost. The present result does not justify a claim
that 214, or any other small number, is sufficient.

The work remains a frozen, finite-basis electronic H6 test. No many-body
sector or new wavefunction was enumerated in discovery or acceptance. The
reference upper's earlier determinant-amplitude discovery is separate. No
general compression theorem, experimental accuracy or physical preparation
capability follows from these intervals.

## Verification and reproduction

**Six focused tests passed.** They check inherited direction and Hamiltonian
bindings, exact target boundaries, rejection of late or weaker results,
complete reconstruction of every full-family generator, and rejection of
the previous diagonal-spin dual by the additional spin-changing PSD gate.
The two selection-boundary tests also passed under `python -S`.

The inherited lower and **all twenty new energy certificates**—nineteen
continuation candidates and one ablation—passed fresh exact replay against
the rational reference upper. No NumPy, SciPy, CVXPY or PySCF was loaded on
that accepting path. All **244 files** in the preceding manifest remain
unchanged.

The [protocol](/Users/aidenlippert/Documents/Spectra/research/spin_enrichment_20260913/PROTOCOL.md),
[proof](/Users/aidenlippert/Documents/Spectra/research/spin_enrichment_20260913/PROOF.md),
and [separate dual protocol](/Users/aidenlippert/Documents/Spectra/research/spin_enrichment_20260913/DUAL_PROTOCOL.md)
record the scope and gates. The
[cost ledger](/Users/aidenlippert/Documents/Spectra/results/spin_enrichment_20260913/cost_ledger.json)
and [fresh replays](/Users/aidenlippert/Documents/Spectra/results/spin_enrichment_20260913/fresh_replay.json)
contain the measured results; proposals, warnings, refusals and hashes are
preserved in
[results/spin_enrichment_20260913](/Users/aidenlippert/Documents/Spectra/results/spin_enrichment_20260913).

From the workspace root, with the configured Python environment:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest \
  research.spin_enrichment_20260913.test_campaign \
  research.spin_enrichment_20260913.test_dual -v
python -S -m research.spin_enrichment_20260913.replay \
  --campaign results/spin_enrichment_20260913/campaign \
  --out /tmp/spin-enrichment-replay.json
```

`campaign.py --out <new-directory>` repeats continuation and ablation with
the frozen limits. `full_dual.py` repeats its separately bounded diagnostic.
Existing discovery output directories are refused. No accepted full-family
witness or physical-error audit is claimed by this pass.
