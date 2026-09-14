# V8 reusable-operation method-chain contract

This contract sharpens the V7 requirement that an acquired procedure `m1`
must causally lower the cost of acquiring a distinct useful procedure `m2`.
It is a test contract, not evidence that the gate has passed. The checker,
raw-data interface, primitive set, budgets, and holdouts must be frozen before
search begins.

## Reusable-operation cost criterion

Let `D` be the same raw public datum, `C` the same immutable certificate
checker, and `P` the same typed public primitive set. For a search policy `S`
define its full acquisition cost on instance `i` as

```text
K(S,i) = K_search + K_observe + K_construct + K_check
         + K_failed/restarted + K_validation,
```

where every primitive call, exact-arithmetic operation, candidate evaluation,
rejection, restart, and checker invocation is charged. Let `Q(S,i)` be a
predeclared quality vector (validity, interval width or residual target, and
refusal status). A reusable operation `m1` qualifies as a method-acquisition
cause for `m2` only if, on prospective paired `m2` instances sampled from
families disjoint from `m1` development, there is a frozen quality constraint
`Q >= q*` and a predeclared aggregate showing a positive downstream effect.
This is an operational acceptance criterion, not a mathematical theorem: its
truth must be established by the preregistered experiment.

```text
NB_H = E[sum_j K(scratch m2_j)
        - sum_j K(m1-enabled m2_j) - K_acquire(m1)]
```

with the same admissibility and quality constraint. `NB_H` is an amortized
net-benefit criterion over a declared horizon `H` of downstream uses. The
separate causal `m2` effect is the paired downstream cost (or valid success at
equal downstream budget) conditional on frozen `m1` availability. Acquisition
may therefore be amortized across preregistered uses; it need not be repaid in
each single `m2` event. A smaller cost on only accepted candidates, or a
faster replay after an uncharged failed search, does not satisfy the criterion.

A positive net benefit is about a complete causal pipeline. A changed candidate
ordering qualifies only when it reduces complete charged work or increases
valid success at equal complete charged work. A useful operation may be
instance-local, but its rule must be reconstructed from `D` and declared
state, not selected from a cache of finished answers.

## Three categories that must be reported separately

**Acquired scientific calculation algorithm.** The search discovers a
reusable, executable rule for constructing or certifying a calculation from
public primitives: for example, a residual-budget allocation rule, a basis
construction with a certified truncation/refusal condition, or a graph-based
partition procedure whose output is checked independently. It is an algorithm
because the rule changes the sequence or representation of valid calculations
and transfers to new families. It must pass the cost criterion and show that the
same rule, rather than a remembered answer, changes `m2` acquisition.

**Ordinary compiler/numerical implementation optimization.** Faster sparse
storage, common-subexpression elimination, vectorization, memoization, better
exact-arithmetic kernels, or a numerically equivalent implementation belongs
here. It can legitimately qualify under the user's allowance for cheaper
evaluation versus better search only when the optimization is itself acquired
from the permitted raw interface, is executable on unseen instances, has its
construction, collision checks, conversion, and validation charged, and
causally reduces the *total* valid `m2` acquisition cost. The optimized
implementation must use exactly the same mathematical primitive/evidence and
the same checker semantics; it cannot quietly use lower precision, extra
observations, a weaker bound, or a precomputed answer. If compiling a
representation lets the second search evaluate more admissible candidates for
the same budget and this produces a held-out improvement, that is a legitimate
cheaper-evaluation method claim. If it merely runs a known supplied method
faster, report it as an implementation optimization and do not call the known
method learned.

**Coefficient tuning.** Choosing thresholds, tolerances, weights, ordering
constants, or a finite hyperparameter from training records is not by itself
algorithm acquisition. It may be a parameter of an already fixed method, but
it cannot satisfy the method gate unless the parameterized rule has a distinct
reusable operation and the full tuned-search and validation cost is charged.
Tuning the checker slack, precision, or refusal criterion to rescue failures
is inadmissible.

The distinction is functional rather than rhetorical: freeze the operation's
typed program and replay it on unseen `D`; then ablate the operation while
holding the representation and all numeric constants fixed. If the claimed
effect disappears only when the operation disappears, and not merely when a
coefficient changes, it is evidence for an algorithmic or implementation
operation. It still requires the cost and causal gates.

## Test contract and anti-leakage requirements

Every arm receives byte-identical raw `D`, the same public primitive/evidence
and precision, the same instance budget, and the same immutable checker `C`.
The checker validates the submitted witness and recomputes all claimed
relations; it never trusts method IDs, search metadata, or cached summaries.
No arm receives a finished partition, optimum, simulator label, answer menu,
hidden target certificate, or evaluator method ID. A residual vector may be
exposed only if it is recomputed from the common generator for every arm.

The experiment must include at least scratch, `m1`-enabled, reconstruction,
restored, and irrelevant-state controls. Reconstruction regenerates the
claimed public summary from `D` with its own charged cost, detecting a cache or
serialization artifact. The restored arm reruns the same disabled acquisition
setup with the byte-identical frozen `m1` procedure and state restored,
establishing reversibility of the claimed cause. The irrelevant arm adds an
executable, equal-size and equal-loading-cost state that performs no useful
operation on the workload and supplies no unavailable evidence; it is not
discovered from unspecified nuisance coordinates. All arms must preserve
failure and refusal paths.

Canonicalize alpha-equivalent programs and no-op reorderings before counting
search. Retain hashes of `D`, the complete trace, submitted certificate,
primitive calls, exact checker receipt, and cost ledger. Report zero-success
instances and negative effects. A control that obtains a known method's output
from the assistant or evaluator is a diagnostic implementation baseline, never
evidence that the method was acquired.

## Conditional candidate two-generation experiment

If headroom and expressivity are first demonstrated, test a candidate
representation algorithm `m1` that constructs a certified
block basis for a Hamiltonian instance. Its typed program may inspect only
observed coefficient magnitudes, local commutation predicates, and residuals;
it greedily grows blocks, orthogonalizes with exact rational operations, and
returns either a basis plus an independently checked omitted-norm budget or
`refuse`. The grammar contains no optimal partition token and no answer list.
Freeze `m1` on development families, including its maximum block size and
primitive budget, after the headroom gate against the frozen coordinate,
residual-adaptive, Arnoldi, and grouping baselines.

Run the second `m2` search on a separate development set to freeze its grammar,
budgets, endpoint, and controls. Then create final `m2` families with unseen coupling patterns, widths, coefficient
scales, and seeds. The distinct second procedure is a residual-budget
allocator: given the same raw `D`, it chooses which remaining coordinate or
Krylov frontier to evaluate next and stops only when `C` verifies the target
certificate. In the enabled arm, `m1`'s frozen block basis is the only
intervention available to the allocator; all basis construction, expansion,
collision checks, and conversion costs are charged. Scratch runs the same
allocator grammar without that basis. Reconstruction rebuilds the basis from
the public summary using charged primitives; restored re-enables the frozen
basis; irrelevant-state supplies matched executable loading overhead with no
useful operator on the workload.

The final evaluation uses common held-out sets `S0`, `S1`, and `S2`: `S0`
measures scratch downstream acquisition, `S1` measures frozen `m1`-enabled
downstream acquisition, and `S2` is the untouched family-level replication.
The paired `m2` gain uses the same instances and checker receipts in `S0` and
`S1`; the horizon net benefit adds `K_acquire(m1)` once across declared uses
and includes both incremental downstream costs and gains. The primary endpoint
is total charged work to obtain a valid certificate with the same predeclared
width/residual target. A secondary endpoint is valid-certificate rate at a
fixed total budget. Analyze effects by independent family. A benefit
that vanishes after basis construction is charged is a null. A benefit that
survives only because the basis changes the checker, precision, admissible
slack, or supplied evidence is invalid. A benefit that survives equal checker
semantics and full accounting is evidence that `m1` acquired a reusable
representation operation which causally made the distinct `m2` search cheaper.
