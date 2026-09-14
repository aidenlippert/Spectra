# Complete-cost test for selective quadratic norm certificates

The diagnostic found that operator squaring can reduce Taylor order on two
existing development cases. This does not establish efficiency. Test a supplied
policy that tries quadratic norms only when the best ordinary bound lies in
`(epsilon, 3 epsilon/2]` and the residual has at most 40 terms. Those thresholds
are development choices, frozen for this test. No held-out outcomes are used.

Compare identical adaptive code with quadratic eligibility capped at one term
(the ablated baseline, with no quadratic trials on this workload) and at 40
terms. Both use the same extended independent checker; the original checker
and prior results remain unchanged. Also retain the earlier adaptive baseline
for context. Use seven paired timing repeats, randomized arm order from seed
271828, on all 36 fixed development problems. Complete elapsed time includes
candidate generation, certificate construction and independent checking.
Failures remain failures and their cost is retained. Report paired effects,
polynomial order, quadratic trial count and acceptance rate. No autonomous
method acquisition is authorized by a norm improvement alone.

For real Pauli coefficients, square exactly by adding the diagonal squared
coefficients to identity and the commuting unordered cross pairs. Anticommuting
cross pairs vanish. If a verified group witness bounds `||R^2|| <= b`, an exact
rational `u >= 0` with `u^2 >= b` bounds `||R|| <= u`. The independent checker
recomputes the square and validates every group. Worst-case pair work is
quadratic in residual support and may outweigh the saved degree.
