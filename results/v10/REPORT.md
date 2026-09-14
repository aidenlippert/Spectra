# V10: exact reference expansion rejected at the complete-cost gate

Compounding scientific capability remains **unproved**. V10 implemented and
checked a conventional solvable-reference calculation, but it does not supply
headroom for method acquisition. No learner, m1/m2 experiment, or reserved
evaluation was run. The next attack concerns an exact matrix implementation of
the existing calculation, described [here](../../research/v10/NEXT_MATRIX_PREFLIGHT.md).

## Fixed workload and outcome

All 36 development inputs retain the same Hamiltonians, central-Z observables,
damping values 0/.2/2, horizons .2/.5 and total operator-norm tolerance .001.
Each arm returns rational Pauli coefficients at the requested endpoint.
Construction, failed attempts, checking and endpoint evaluation are timed.

| Supplied conventional method | Certified | Complete seconds, including refusals |
| --- | ---: | ---: |
| Fraction-free Taylor | 26/36 | 2.345 |
| Strengthened BFS | 27/36 | 4.501 |
| Solvable-reference expansion | 21/36 | 82.634 |

The reference arm loses on all 21 common successes. Its median cost ratio is
35.10 relative to Taylor and 27.68 relative to BFS; even its best relative
case costs 10.67 times Taylor and 5.34 times BFS. It adds no certified case.
This is one diagnostic pass, not a general timing theorem. The large losses
and reduced coverage are sufficient to close this implementation without a
positive timing-confidence campaign.

Fifteen reference refusals comprise six live-mode support failures, six cached
column failures and three total-entry failures. Refusal work costs 44.19 seconds
and remains in the total. A budget refusal is not a physical impossibility
certificate.

## What was actually constructed and checked

Split G=G0+B, where G0 contains existing onsite-Z fields and the same local
depolarization. All other nontrivial Hamiltonian terms remain in B=i[V,.].
The exact I/Z/P/M eigenbasis permits rational complex exponential-polynomial
corrections Dk. Coincident frequencies produce polynomial factors through exact
integration; no numerical near-equality or singular division is used.

The independent checker verifies D0, every correction's differential identity
`Dk'=G0 Dk+B D(k-1)`, its zero initial condition, and conjugate Hermiticity.
It derives K=2 sum|V_P| and verifies the model-conditional factorial bound
`||O||_1 (KT)^(m+1)/(m+1)!`. Here `||O||_1` denotes the sum of absolute input
Pauli coefficients, an upper bound on operator norm. Both physical semigroups
are contractive on Hermitian operator norm, so this estimate requires no
additional exp(KT) factor. These are supplied mathematical tools, not new
discoveries by the system.

Endpoint evaluation shares equal exponent calculations, aggregates exact
polynomial coefficients first, and uses rigorous rational Taylor disks with
scaling, squaring and charged dyadic rounding. The allocated truncation error
is at most .000999 and readout error at most .000001. Exact conversion to
Pauli coefficients and taking real parts implements Hermitian symmetrization;
the [interface proof](../../research/v10/readout.md) shows that the operator-norm
error cannot increase.

Among accepted reference calculations, construction costs 19.33 seconds,
checking 16.80 seconds and endpoint evaluation 2.31 seconds. There are 170,444
retained layer entries across those cases. The symbolic construction and its
checking dominate; simply improving exponential evaluation would not remove
the observed loss.

## Verification and corrections

- All **74 exported certificates and endpoint outputs** replayed: 26 Taylor,
  27 BFS and 21 reference. The archive is approximately 14 MB.
- **180 tests passed**, including exact basis round trips, dense commutator
  checks, recurrence and initial-condition tampering, repeated frequencies,
  independent rational-series enclosures, resource refusals, and a two-qubit
  physical endpoint comparison. Tests are engineering verification, not
  evidence of compounding capability.
- Forty historical source/data/result hashes still match. V7/V8/V9 archives,
  V8/V9 sources and the original V7 checker remain unchanged. V4 calibration
  remains outside V5's dependency ledger.
- Review rejected an incorrect agent exponential remainder explanation and
  inadequate tolerance allocation. The readout was rewritten and checked.
  An agent dense test also compared B with full G's Hamiltonian; it was fixed
  to compare the actual interaction. Invalid drafts are retained under
  `research/v10/invalidated/` and are not used as evidence.

[Frozen protocol](../../research/v10/HEADROOM_PROTOCOL.md) ·
[Complete measurements](headroom.json) ·
[Replay receipt](verification_receipt.json) ·
[Exact archive](certificate_archive.json) ·
[Checker proof review](../../research/v10/root_checker_audit.md).

## Decision

Retire this exact reference-expansion schedule for the tested workload. Its
mathematical validity does not imply algorithmic advantage. The broader goal
remains active: a feasible construction improvement must precede acquisition,
and any future m1 must causally reduce m2 acquisition cost before two untouched
incremental gains and acquisition-inclusive benefit can be claimed.
