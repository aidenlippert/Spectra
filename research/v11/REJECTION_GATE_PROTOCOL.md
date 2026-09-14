# V11 conventional rejection-gate diagnostic

The matrix comparison failed. Its profile identifies repeated norm-proposal
work as a larger remaining cost than the matrix action itself. For Hermitian
R=sum r_P P, the exact inequalities

`||R||infinity >= sqrt(sum_P r_P²) >= max_P |r_P|`

give necessary conditions for a valid residual norm upper bound to pass.
The first follows from normalized Frobenius norm and Pauli orthogonality.
At Taylor order k, if the weighted lower bound exceeds the fixed tolerance,
skip the ordinary l1/firstfit/weighted proposals. Strict equality proceeds to
the original proposals. This preserves the first accepted polynomial and norm
witness, conditional on the same arithmetic resource bounds.

Three predefined schedules are compared against frozen V8: max-coefficient
only, squared normalized Frobenius only, and max followed by Frobenius. Every
guard operates on the existing integer derivative numerator before creating
Fractions for futile norm proposals. Early exits and inspected/squared entries
are recorded. The original dynamics, degree24 cap,512 Pauli/cache budgets,
grouping algorithms and V7 checker are unchanged.

First check the existing n3XXZ gamma2 T.5 case for actual complete-cost benefit.
If a guard clearly improves it, use the same36 development inputs in five
randomized repetitions, seed111. Include all failed attempts and fresh checking
plus endpoint evaluation. Record source hashes and preserve every failure.
For every accepted input, require exact equality of polynomial, norm witness,
bound and endpoint across arms. Use paired timings and total costs including
refusals; do not select a favourable physical subset or change damping.

These are supplied conventional inequalities and schedules. A speedup strengthens
the baseline and is not an acquired m1 or a compounding result. A future learner
must have these operations available in its matched conventional comparisons.
