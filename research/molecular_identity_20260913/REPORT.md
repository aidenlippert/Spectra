# Molecular local-operator pass — 13 September 2026

The pass found a useful compression mechanism on H4, but the same selection
rule did not deliver accurate H6 bounds. Coupling short local cubic operators
to the existing global operator blocks narrowed H4's exact interval from
3.563811 to **0.537323 millihartree**. H6 improved only from
12.208352 to **11.990119 millihartree**. The target is 1.6.

This advances a specific certificate-design question within the broader mission
of affordable, trustworthy prediction and control of matter. It establishes no
new material, experimental capability, or general solution of the many-body
problem.

## Outcome and controls

All widths below are exact rational upper-minus-lower intervals, shown in rounded
millihartree. Both ends refer to the same frozen rational Hamiltonian. H4 has
eight spin orbitals/four electrons; H6 has twelve/six. They are the existing
STO-3G hydrogen-chain controls.

| Method | H4 width, mHa | H6 width, mHa | H4 / H6 accuracy target |
|---|---:|---:|---|
| Previous creator-channel certificate | 3.563811 | 12.208352 | Fail / Fail |
| All three-orbital local blocks | 3.642096 | 12.328762 | Fail / Fail |
| Separate local blocks plus two rounds of selected four-orbital blocks, best exported result | 3.091607 | 12.215911 | Fail / Fail |
| Short cubic directions coupled into existing blocks, best exported result | **0.537323** | **11.990119** | **Pass / Fail** |

H4 first crossed the target after one enrichment round at 1.497423 mHa and
3,704 Gram entries. Its final result retained 36 added generators including
adjoints, each supported on at most four spin orbitals. The added H4 cubic
polynomials contain 5–12 monomials. The H6 campaign retained 48 generators;
its best result occurred after round two. Round three gave 11.994490 mHa,
slightly weaker after numerical solving/export. The best previous certificate
remains available whenever a new exported result is weaker.

A direct H4 control reused the **same 36 learned generators** and placed them in
separate PSD blocks, removing their cross terms with the base. Its interval was
**3.429122 mHa**, versus **0.537323 mHa** when coupled. This supports the practical
importance of those cross terms for this dictionary. These are comparisons of
accepted exported certificates, not proofs of the exact optima of the two cones.
The ablation inherits the 10.55-second dictionary-discovery cost and adds a
3.67-second solve/export; that inherited cost cannot be omitted.

## Exact structural finding

The previously accepted restricted H4 dual assigns the occupation projector
`n₂ n₄ n₅` the exact value

`-87413147645264941 / 7000000000000000000`.

CAR gives `(a₂ a₄ a₅)†(a₂ a₄ a₅) = n₂ n₄ n₅ >= 0`. The negative value therefore
proves that this restricted dual violates a valid physical positivity constraint.
The eight three-orbital occupation projectors sum to identity, independently of
particle number. Their validity transfers to any distinct orbital labels without
assuming a fixed number of electrons inside a cluster.

The initial scan inspected 56 three-orbital and 70 four-orbital subsets, for
charge minus one and minus three. Its discovery scan took 1.08 seconds, plus
2.78 seconds to check the old dual. A separate standard-library replay accepted
32 stored rational negative directions and all eight occupation identities.
The earlier exact-dual construction itself took 7.28 seconds in its original
receipt; it is reused diagnostic evidence, not a free new discovery step.

Eliminating one nonphysical dual does not establish a better optimum: another
feasible dual may replace it. The all-three-orbital solve demonstrates why the
structural diagnostic must be followed by actual certificate construction.

[PROOF.md](/Users/aidenlippert/Documents/Spectra/research/molecular_identity_20260913/PROOF.md) gives the operator identities, coupling argument, and exact
lower/upper acceptance conditions.

## Fixed rule tested on H6

The coupled method starts afresh from the Hamiltonian and the creator-channel
cone. In each of three rounds it examines all four-orbital subsets, partitions
linear/cubic operators by exact charge and Hamiltonian symmetries, and ranks
negative numerical moment directions. It selects at most eight, rounds their
coefficients to rationals with denominator one million, removes the already
represented linear part, and adds the cubic parts and their exact adjoints.

An existing block with the same charge/symmetry receives the new generators,
allowing global cross terms. A previously absent channel receives a new PSD
block. H6 exposed a missing implementation case for that latter operation before
round-two construction. The exception was preserved, the case was covered by a
focused test, and the saved run resumed. H4 required no such new channels. The
selection policy, thresholds, round count, and solver budgets were unchanged.

The exact H4 witness informed the research design. Neither coupled discovery
run reads its coefficients, old certificate factors, nor the reference state.
H4 and H6 rebuild their proposals from their own Hamiltonian and numerical duals.
This is a transferred **selection rule**, not a claim that one set of learned
H4 coefficients predicts H6 accurately.

## Complete resource accounting

| Coupled campaign | H4 | H6 |
|---|---:|---:|
| Discovery including seed, all rounds, pricing, maps, solve, export/check | 10.55 s | approximately 270.74 s |
| Map/construction time summed over solves | 1.44 s | 27.22 s |
| Solve time summed over solves | 7.49 s | 209.11 s |
| Export plus embedded exact lower verification | 1.11 s | 29.58 s |
| Local-direction pricing | 0.37 s | 3.93 s |
| Pricing word-pair products, all rounds | 42,792 | 330,264 |
| Additional independent interval replay, all solved stages | 1.62 s | 22.14 s |
| Gram scalar entries at best stage | 4,484 | 17,868 |
| Coefficient-map nonzeros at best stage | 83,410 | 1,363,308 |
| Best certificate bytes | 179,855 | 1,044,595 |
| Best certificate factor rows / nonzeros | 239 / 7,624 | 501 / 48,536 |
| Factor denominator bits / maximum numerator bits | 68 / 68 | 68 / 67 |

The H6 total includes a 121.08-second pre-failure wall-time estimate recovered
from timestamps. Time spent diagnosing the code failure is not included as
computation. Some independent runs overlapped, so these are observations from
one campaign rather than controlled speed comparisons. Solver limits are soft
60-second requests per solve; construction, canonicalization, and export add
cost, and observed solve wall time can exceed that request.

The weaker separate-local-block campaigns consumed 9.52 seconds on H4 including
its reused seed, and 257.20 seconds on H6. Their costs and unsuccessful intervals
are retained. The cost ledger is [summary.json](/Users/aidenlippert/Documents/Spectra/results/molecular_identity_20260913/summary.json).

For context, the previous larger chosen cubic dictionary used 17,776 H4 and
217,268 H6 Gram entries and obtained much tighter widths: 0.000153 and 0.054994
mHa. Its recorded source wall times were 7.26 and 259.31 seconds, with different
jobs/environments and additional proof/reference costs. The new dictionaries
use fewer PSD entries, but this pass demonstrates neither a speed advantage nor
smaller certificate files: the older files were 164,227 and 970,728 bytes.

The existing upper witnesses contain 20 H4 and 200 H6 nonzero determinants.
Their original FCI discovery receipts report 36- and 400-dimensional spin
sectors and 0.015/0.233 seconds at these tiny sizes. Upper replay performs
3,680 and 183,600 word/state checks. That exponential reference discovery is
kept outside lower-bound discovery and explicitly charged as reference work.
For these small controls FCI is far cheaper; this experiment probes certificate
structure rather than practical superiority over exact diagonalization.

## What is established and what remains open

All **15** exported intervals were independently replayed with rational
arithmetic. **18 focused tests passed**, including unchanged solver controls,
CAR cross terms spanning more than four orbitals, occupation identities,
corrupted witness rejection, Hamiltonian binding, export degree-six retention,
and the newly discovered symmetry-channel case. Four exact-only tests also ran
under `python -S`. Numerical solves reported `optimal_inaccurate`; their status
is never the accepting criterion. The exported factors and exact residual give
the reported lower bounds.

At the best H6 stage the exact residual penalty is only
0.042267 mHa, around
0.35% of its remaining interval.
Eliminating that penalty would still leave it far above the target. Numerical
optimization of the chosen cone may still matter: no exact H6 dual obstruction
for the enlarged dictionary has been proved. The outcome is a failed finite
transfer experiment, not an impossibility theorem.

Selecting local positivity constraints is established RDM/SOS methodology;
see [Li–Lu's cluster-selection work](https://arxiv.org/abs/2305.18571) and the
[primary-source audit](/Users/aidenlippert/Documents/Spectra/research/molecular_identity_20260913/SOURCES.md). Generic cubic squares use degree-six moments;
they are not automatically the reduced T1/T2 constructions. No new physical law
or general compression theorem is claimed.

**Next course:** retain the coupling mechanism and focus on H6 under a fixed
end-to-end budget. Test whether ranking small sets of directions by actual
energy-bound improvement per construction/solve/replay cost outperforms ranking
only by the most negative local moment. Use an exact dual obstruction if claiming
that the current span itself is insufficient. Larger four-to-six-orbital motifs
remain candidates, but adding them is justified only by measured bound gain.
The present selection rule reduced H4's interval by about
84.9% and H6's by only
1.8%.

Reproduction commands are in [REPRODUCE.md](/Users/aidenlippert/Documents/Spectra/research/molecular_identity_20260913/REPRODUCE.md). Existing fixtures,
upper witnesses, and previous results were preserved. The reusable solver change
is a validated optional additional-block/merge interface; its original source
and diff are saved with the campaign artifacts.
