# V7 method-learning gate review: many-body certificates

This is a prospective adversarial review, not a result. The target is a
two-generation transfer: acquire a research procedure `m1` from raw public
problem data, then let `m1` causally change acquisition of a second procedure
`m2`. The claim is allowed only if `m1` lowers the total cost of discovering
`m2` on prospective holdouts, with equal certificate validity. A good
certificate produced by a supplied method is not acquired method transfer.
A discovered faster implementation can qualify as an acquired computational
method if the causal and cost gates below are met; it is not automatically
improved hypothesis selection.

## Freeze the first gate before method search

The fixed evaluator should generate many-body Hamiltonian instances and a
public primitive representation from the supplied common physical prior. The
current `experiments/certificates.py` is an energy-certificate checker, not a
dynamics or residual checker. A new dynamics/residual checker may be used only
after its source, hash, exact arithmetic rules, and refusal tests are frozen
and independently verified. Every candidate gets the same raw datum `D` and
declared compute budget. The evaluator owns held-out validation instances and
the common-prior generator; candidates do not receive a finished partition,
an answer menu, simulator labels, or a hidden target certificate.

Before any training, freeze the physical depolarization rates `gamma = 0, .2,
2`, the instance-generation seeds, and the candidate-evaluation budget.
Gamma must enter the supplied physical observation/noise model identically for
every arm; it is not certificate slack and cannot be tuned after seeing a
method's failures.
For each raw instance retain an immutable hash of `D`, the candidate trace,
the submitted certificate, and the checker receipt. Re-running the checker on
the same receipt must be independent of the search procedure that created it.

The baseline gate is a strong, fixed-rational residual checker with a frozen
headroom report. Compare at least the following supplied strategies at all
three gamma values: coordinate Krylov expansion, residual-adaptive expansion,
Arnoldi-style orthogonalization, and group-norm certificate grouping. These
are baselines and diagnostic arms, not training labels. A candidate cannot be
declared useful because it beats a weak coordinate order or because its output
is accepted only at a larger unsupported slack.

Report certificate validity, absolute lower/upper interval width, residual
bound, refusal rate, number of primitive observations, search operations,
checker operations, and evaluator cost separately. Charge failed candidates,
restarts, rejected certificates, and validation evaluations. Also report
cost per candidate evaluation and cost per retained candidate; averaging over
only successful candidates would reward early rejection or hidden retries.

## A minimally expressive search language

The search object should be a small typed program over public primitives, not a
library of named finished methods. A useful core has:

* selectors over observed coordinates, residual coordinates, and declared
  commutation/anticommutation predicates;
* bounded `map`, `fold`, `partition`, `orthogonalize`, `truncate`, and
  `refuse` operators;
* exact rational arithmetic plus explicitly bounded norm and square-root
  operations;
* a state record containing only prior observations, residuals, accepted
  bounds, and charged work; and
* a final constructor for a certificate witness that the immutable checker
  validates.

Programs have a fixed maximum AST size, depth, and primitive-call budget. The
grammar may express coordinate order, residual thresholds, Krylov updates,
Arnoldi reorthogonalization, and group-norm aggregation, but it must not
contain a token meaning “use the optimal partition”, a hidden certificate,
the evaluator's method ID, or a supplied answer list. Constants are either
publicly fixed or learned from training records with a separate validation
split. Alpha-equivalent programs and no-op reorderings should be canonicalized
before counting candidates, so syntactic duplication cannot manufacture a
large search history.

The evaluator should expose outcomes, not oracle explanations: accepted or
refused status, checked rational inequalities, and charged work. It may expose
the residual vector because residual-driven research is the object of study,
including omitted terms recomputed from the common supplied generator; it
must not expose held-out optima or privileged model information. A
program that cannot establish its own error budget must return `refuse`.

## Causal two-generation design

Use disjoint prospective families for development, `m1` acquisition,
`m2` acquisition, and final evaluation. A generation consists of multiple
independent Hamiltonian families, not repeated seeds of one family. The
training set may teach `m1`; the second-generation set is created after `m1`
is frozen and contains new coupling patterns, widths, coefficient scales, and
unseen seeds. Hold out an entire family type for the final test, so retrieval
of an instance signature cannot masquerade as transfer.

The causal intervention is the only difference between the `m2` arms:

1. **scratch:** search `m2` from raw `D` under the frozen grammar and budget;
2. **m1-enabled:** run the frozen learned `m1` on the same raw `D`, charge its
   full cost, and permit its returned procedure/state to change the next
   observations, ordering, compression, or rejection policy;
3. **reconstruction:** disable the learned state and reconstruct the same
   public summary from `D` using the declared primitive, to detect a label or
   serialization artifact;
4. **restored:** restore the exact frozen `m1` procedure and verify restoration
   of its downstream effect;
5. **irrelevant-control:** provide an equally sized state computed from
   nuisance coordinates that cannot affect the certificate family.

The m1-enabled arm must improve a predeclared downstream quantity (for
example, total primitive work to reach a fixed certificate width and validity
target, or the probability of obtaining a valid certificate within a fixed
budget). A change in candidate ordering counts only if the full search and
verification costs fall. If it merely finds the same candidate earlier while
the charged total is unchanged, the result is null. If `m1` changes the
representation, the evaluator must show that this change is what causes the
`m2` improvement by replaying the same raw `D` under the disabled and restored
interventions.

## Critical gates

The following gates are conjunctive.

1. **Checker gate.** Corrupted coefficients, duplicate/missing terms,
   invalid clique relations, false square-root intervals, and inconsistent
   omitted norms are rejected. The checker must never trust search metadata.
2. **Primitive expressivity gate.** Every claimed improvement is executable
   using the public typed primitives and equal raw `D`; no method receives
   extra observations, stronger precision, or an evaluator-only feature.
3. **Baseline headroom gate.** At each fixed depolarization rate, at least one candidate
   family must beat the strongest frozen baseline on a predeclared paired
   cost/quality measure. If none does, stop method acquisition for that family
   and develop a bounded new candidate family; do not tune the search toward
   a convenient null.
4. **First-generation gate.** `m1` must be frozen using only its development
   split and must have a reproducible rule, not a cached table of instances.
   It must survive family and seed holdouts and have no greater unsupported
   claim or refusal defect than the strongest baseline.
5. **Causal transfer gate.** On untouched `m2` families, m1-enabled must beat
   scratch after charging m1 acquisition and all downstream evaluations. The
   restored control must restore the frozen `m1` state and its claimed causal
   effect; it is not a placebo or randomized state. The irrelevant control must
   be unable to reproduce that effect. A reconstruction arm may recover an
   equivalent benefit if its reconstruction cost is charged; it is a control
   for serialization/cache artifacts, not a required null.
6. **Prospective holdout gate.** The final family is never used for grammar,
   hyperparameter, stopping, or gamma selection. Its exact checker receipts,
   raw-data hashes, and candidate traces are retained for replay.

No gate is passed by a single successful certificate. Use multiple independent
families and seeds as a prospective replication design, with the exact number,
paired effect threshold, and uncertainty procedure frozen before evaluation.
Three families and two seeds per family are a useful minimum for detecting
obvious brittleness, not a statistical confidence guarantee. Report the effect
distribution over the true sampling units, including zero and negative cases;
candidate counts are not independent observations.

## Candidate reusable procedures

The most credible reusable objects are meta-procedures that learn *when and
where to spend work*, while the checker remains fixed. They can be expressed
in the typed language above without naming a finished certificate.

* A residual-budget allocator can estimate which coordinate block or Krylov
  frontier has the greatest certified reduction per primitive operation, then
  stop or refuse when the reduction is below its measured cost. Training can
  learn the allocation rule; the held-out checker still validates every
  bound.
* A commutation-graph ordering rule can learn a permutation or partition
  heuristic from graph statistics (degree, degeneracy, local coefficient
  scale, and residual history). It must be compared with greedy grouping and
  exact small-graph optima only as evaluator diagnostics, never supplied as
  labels to the learner.
* A truncation-and-refusal rule can learn a certified threshold for dropping
  weak or rapidly decaying residual coordinates under each physical gamma.
  The dropped mass must enter the independently checked omitted-norm budget;
  a predictive score alone is not evidence of a valid certificate.
* A cross-instance compression rule can detect repeated local operator
  neighborhoods and choose a reusable basis or block representation. Its
  benefit counts only if basis construction, collision checks, expansion, and
  checker work are cheaper in full than rebuilding each instance.

These procedures are plausible because they can reduce complete acquisition
and verification work across changing instances without encoding which
partition or certificate wins. They remain hypotheses: each must first show
headroom against all four frozen baselines, then show that its frozen rule
changes the second-generation search on a family held out by coupling pattern
and width. A learned order that leaves the number of checked terms and exact
arithmetic unchanged has no net method benefit.

## Adversarial failure modes

The most likely false positives are leakage of the held-out optimum through
residual normalization, a grammar primitive that encodes one baseline's
answer, unequal rational precision, uncharged failed searches, and caching
that benefits only the enabled arm. Also test scale changes, permuted
coordinate names, duplicated terms, near-degenerate residuals, dense
anticommutation graphs, omitted-term misspecification, and cases where no
admitted certificate reaches the target. A valid method must abstain on the
last case rather than convert an unverified residual into a certificate.

If all candidates fail, preserve the null, identify whether the obstruction is
expressivity, headroom, identifiability, or evaluator cost, and redirect only
to a bounded new candidate family with a new prospective gate. Do not mark the
goal complete until the complete frozen chain `m1 -> m2`, the causal controls,
and net benefit after all acquisition and evaluation costs have passed.
