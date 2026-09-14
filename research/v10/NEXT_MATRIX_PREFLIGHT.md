# Next bounded preflight: exact matrix action and fast Pauli conversion

V10's reference expansion has failed the development cost gate. Do not tune
its mode budget or reopen it with favorable damping. A distinct pending cost
question concerns the representation used to execute the existing Taylor
calculation. This is a conventional baseline question, not acquired m1.

## Physical matrix action

For d=2^n, retain O as a d-by-d complex matrix, using separate exact integer
real/imaginary arrays after clearing denominators. Apply

`G(O)=i[H,O]+gamma/4 sum_(i,a=X,Y,Z)(sigma_ia O sigma_ia-O)`.

This gives exactly `-gamma weight(P) P` on each Pauli P. Multiplication or
conjugation by each Pauli term is a signed permutation with phases among
{1,-1,i,-i}. Applying L Hamiltonian terms and 3n damping conjugations costs
O((L+n)d²) entries per action. Alternatively a dense H costs O(d³+n d²).
Neither route requires a d²-by-d² Liouville matrix. Storage is O(d²), still
exponential in n; only n3/4/6 are development inputs.

## Exact conversion without trace-by-trace matrix products

For bit masks x,z, let `y=popcount(x&z)` and
`P(x,z)=i^y X^x Z^z`, with each qubit letter I,X,Z,Y. Then

`c(x,z)=trace(P(x,z) A)/d`
`       = i^y/d sum_s (-1)^popcount(z&s) A[s,s xor x]`.

Proof: `P|s>=i^y(-1)^(z·s)|s xor x>`; expand the trace and substitute.
For each x, the inner sum over s for all z is a Walsh-Hadamard transform.
Thus conversion uses O(d² log d) additions/subtractions and exact phase/division
operations. This is established baseline mathematics: see
[Georges, Berntson, Sünderhauf and Ivanov](https://arxiv.org/abs/2408.06206)
and the related [tensorized algorithm of Hantzko, Binkowski and Gupta](https://arxiv.org/abs/2310.13421).
Neither the transform nor this rederivation may be claimed as autonomous discovery.

The coefficient formula was checked against direct traces for every Pauli
coefficient on seeded small integer-complex matrices n1..4 (4+16+64+256=340
equalities, exactly representable binary rationals). This is an arithmetic
sanity check, not a performance or scientific-transfer result; see
`results/v10/pauli_transform_sanity.json`.

## Strongest narrow first implementation

Initially implement the *same exact* fraction-free Taylor polynomials as V8
using matrix action plus fast Pauli conversion, without introducing roundoff.
The original V7 checker and the existing ordinary norm grouping rules can then
be used unchanged. After each integer derivative, convert to Pauli coefficients
with all work charged, construct the same tail witness, and include endpoint
conversion and fresh checking. Verify exact coefficient equality to V8 for
small cases before any timing claim.

Use int64 only when an integer forward bound proves all products, partial sums
and Walsh butterflies fit. For a signed-permutation interaction implementation
with cleared coefficients h_q and g_q=q gamma/4, a safe componentwise action
bound is `(2 sum|h_q|+6n|g_q|) M` for a matrix whose real and imaginary entries
have magnitude at most M. A Walsh row sum is bounded by d M. Check bounds with
Python integers before entering a fixed-width kernel. If either bound fails,
use an explicitly costed arbitrary-integer fallback or refuse; never detect
overflow only after it occurred. Track temporary magnitude as well as final
cancelled magnitude, rational bits, matrix entries, and Pauli support caps.

First paired case: existing n3 XXZ, gamma2, T1/2, epsilon1/1000, central Z.
If the complete method loses, retain the exact conventional result and identify
whether action, conversion or immutable checking dominates before considering
another implementation. A win only licenses development expansion and proper
paired repetitions; it does not open heldout evaluation or establish learning.
