# V8: exact norm covers, stronger arithmetic, and rejected headroom candidates

**The two-generation scientific-method goal remains unachieved.** No m1/m2
learner was launched and no reserved heldout outcomes were generated. V1–V7
remain preserved, including the actual-model TCLab refusal gate and the V6
history-state rejection. V4 calibration remains outside V5's dependency ledger.

This pass used independent agents for cost profiling, coefficient covers,
cache lifetime, algebraic structure, time stepping, and the causal experiment
contract. Root review corrected mathematical, baseline, and accounting errors
before accepting results. Agent completion messages were not treated as proof.

## Complete-cost development replay

[Protocol](../../research/v8/COST_PROTOCOL.md),
[all 720 rows](cost_replay.json), and
[exact constructions](../../research/v8/CONSTRUCTIONS.md).

All arms use the same 36 development physical cases: widths 3/4/6, XXZ/mixed
Hamiltonians, fixed gamma 0/.2/2, horizons .2/.5, tolerance .001, live support
and cache caps 512, degree cap 24. Five repetitions per case randomize arm
order. Timings include fresh construction, every attempted norm computation,
solve/exports, independent checking, and refusal work. Import startup is shared
and outside per-call time. These are development timings, not independent
samples from a heldout distribution.

| Arm | Certified per repetition | Refused | Total seconds across five repetitions |
|---|---:|---:|---:|
| Conventional rational Taylor | 26/36 | 10/36 | 11.488 |
| Conventional integer Taylor | 26/36 | 10/36 | 11.005 |
| Partition Taylor through shared V8 norm checker | 26/36 | 10/36 | 11.847 |
| Overlapping-cover Taylor through the same checker | 26/36 | 10/36 | 11.858 |

Integer Taylor had a median paired speedup of **1.060×** on the 130 successful
paired calls, winning 127 of them. Total cost including refusals fell about
4.2%. Exact polynomials and certificate bounds remain identical. This is a
useful implementation improvement using supplied denominator-clearing math;
it strengthens the baseline and does not qualify as acquired m1.

The overlapping cover reduces Taylor degree 12→11 in the n3 XXZ, gamma2,
T=.5 case. Across its five paired timings it is **2.09–2.25× slower** despite
that smaller polynomial. All other cases retain the same order and feasibility.
Whole-suite timing difference is small because the costly cover runs rarely;
there is no positive headroom result. This candidate is closed in its current
form. Retaining fewer terms or using one fewer order is not the criterion.

## What the norm mathematics did establish

The new exact checker validates overlapping anticommuting coefficient splits.
The five-cycle Pauli example gives a strict bound improvement over **every**
disjoint partition: 5/√2 versus 1+2√2. The proof and coordinate minimizer are
in the constructions document. No priority claim is made.

On 48 actual development residuals the bounded proposer tightened all 48
bounds, but only one additional residual crossed tolerance. Cover construction
and checking cost substantially more than ordinary grouping. Every split and
its exact residual coefficients is saved in [the norm archive](fractional_probe.json).
The mathematical improvement does not survive the full-cost gate.

## Other completed attacks

- **Cache lifetime:** a bounded LRU policy with unchanged live support recovered
  zero cases. It preserves 26 successes and ten refusals. The observed cache
  refusal later becomes a live-support refusal. A larger support cap would be
  a different resource budget.
- **Time stepping:** 172 order attempts over all 36 cases and three segment
  counts produced 72 certified attempts, 64 over-tolerance results and 36
  refusals. Each segment arm succeeds on 24 cases, losing the two width-6
  successes available to adaptive Taylor's finer degree schedule. No new
  feasible case appears. Extra segment costs are paid; earlier residual error
  remains accumulated even when interface jumps are exactly zero.
- **Operator structure:** the reachable commutator graph is bipartite through
  Y-count parity. A linear-charge search missed this nonlinear charge. This
  known real/imaginary structure already removes forbidden sparse edges; no
  additional operation reduction has been established.

These are bounded negative or conventional results, not impossibility theorems.
Reports are in [research/v8](../../research/v8/).

## Review corrections and evidence status

The first cover prototype used nearest floating square roots, produced invalid
cliques, and compared a triangle cover against an incomplete partition baseline.
It was replaced, and its artifacts were preserved under
[invalidated](../../research/v8/invalidated/). The corrected checker verifies
exact coefficient reconstruction and outward squared-root inequalities.

The initial time-step loop skipped physical cases, then a revision overwrote
an accepted baseline witness with an inferior l1 witness. Both errors were
corrected before the final all-case replay. The profile's attribution of a
`norm_witness` call to the independent checker was rejected by call-graph
inspection. No sampled recurrence-checking shortcut was introduced.

The independent replay accepted 104 evolution certificates and 48 norm
certificates, verified 26 exact integer/rational polynomial pairs, and replayed
40 refusal outcomes. All 162 tests pass. Forty historical source/data/result
hashes match and the V7 certificate archive is unchanged. The receipt is in
[verification_receipt.json](verification_receipt.json). Passing tests and exact
replays establish the tested implementation's certificate behavior; neither
is counted as scientific learning or hardware validation.

## Next attack

The remaining target is a cheaper **construction of the whole approximation**,
not another tiny norm-threshold variation. The next bounded preflights examine
non-Taylor polynomial construction and trajectory-wide residual-driven basis
expansion, with the integer baseline and standard spectral/Galerkin mathematics
available to all relevant arms. Each must first demonstrate feasible complete
cost headroom. The method-search and two-generation gates remain unchanged;
[the causal contract](../../research/v8/method_chain_contract.md) now explicitly
separates acquisition cost, amortized benefit, evaluator speed and search quality.
