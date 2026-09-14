# V12 bounded audit: local Pauli-action signature compilation

## Finding

The strongest remaining repeated operation is exact sparse generator action. In
`experiments/v11_guarded_taylor.py`, `apply_integer` requests one column per
currently present Pauli label, scans every Hamiltonian pair on a cache miss,
then performs coefficient multiply-adds and dictionary merges over the column
(`apply_integer`, especially the `commutator_i(gen.h, p)` path). The Frobenius
guard removes impossible norm proposals, but it does not remove these actions;
the report says that all retained polynomial coefficients and generator actions
are still constructed (`results/v11/REPORT.md`, “Mechanism and proof”). V9's
prefix/difference diagnostic only measured reuse of already computed trajectory
prefixes and found no complete-cost win (`results/v9/REPORT.md`, “Incremental
reuse diagnostic”); this proposal is therefore a different representation
operation.

## Candidate reusable operation

Compile an exact **local Pauli-action signature table** from the public
Hamiltonian `h` before the second calculation. For each Pauli word `p`, the
table records the ordered local rewrite events needed by `commutator_i(h,p)`:
the term index, overlap/commutation result, output Pauli word, phase/sign, and
the damping contribution. The compiler is allowed to inspect only `h`, `p`,
the declared Pauli primitive, and the fixed support/cache budgets. It must
canonicalize alpha-equivalent local events and emit either a typed signature
program or `refuse` when the table or program exceeds a fixed budget. Applying
the program still uses exact integer arithmetic and must reproduce the ordinary
column exactly; no norm, checker, or physical parameter may be changed.

The useful two-generation hypothesis is that `m1` constructs this reusable
local rewrite representation from each public instance, while `m2` is a distinct residual or
Taylor search that requests many columns. The enabled arm uses the frozen
signature program; scratch constructs columns directly; reconstruction rebuilds
the signatures from `h` with charged work; restored re-enables the same frozen
program; irrelevant-state loads equal-size executable state with no column
semantics. Every construction, collision check, serialization/conversion,
failed compilation, and exact checker call is charged. The checker remains the
unchanged V7 checker (`experiments/v7_certificate.py`).

This is falsifiable on the current physics: for a disordered ladder, coupling
values may differ while the local Pauli rewrite topology remains compilable,
so the signature must be generated from each instance rather than copied from
finished answers. The observable output is the exact same Taylor polynomial,
residual witness, and endpoint. A mismatch, unsupported topology, or any
reliance on a method ID/cache is a refusal, not a success.

## Cost-based disqualifier

Close the candidate if, on a frozen paired development family, charged
`K_compile + K_apply + K_check` is not lower than direct exact action at the
same valid-certificate target, or if the gain disappears after reconstruction.
For eventual acquisition, also require a positive horizon net benefit under
the V8 contract:

`NB_H = E[sum scratch m2 cost - sum enabled m2 cost - K_acquire(m1)] > 0`.

The candidate is especially likely to fail when each column is used only once,
when signature construction scans nearly all Hamiltonian pairs, or when map
merge work dominates. No isolated column-count reduction or faster replay is
evidence of acquisition. The first bounded test should use only existing
development cases and compare complete charged cost against the V11
Frobenius-guarded Taylor baseline; it must not open the reserved width-5/7,
disordered-ladder central-X/Y evaluation.

## Exact source locations and boundaries

- `experiments/v11_guarded_taylor.py`: integer recurrence, column cache,
  Hamiltonian-pair scan, and coefficient multiply-add accounting.
- `experiments/v7_certificate.py`: authoritative `Generator.apply`, Pauli
  action, and immutable certificate checker.
- `results/v11/REPORT.md`: Frobenius guard's measured cost and the fact that
  generator actions remain charged.
- `results/v9/REPORT.md`: prefix/difference reuse diagnostic and its null
  complete-cost result.
- `research/v8/method_chain_contract.md`: required scratch, enabled,
  reconstruction, restored, irrelevant controls and horizon accounting.
- `research/ACTIVE_STATUS.md`: reserved evaluation and no-acquisition status.

This audit identifies a mechanism worth one bounded compiler-style preflight;
it is not a claim that the operation has been acquired or that it improves the
scientific calculation.
