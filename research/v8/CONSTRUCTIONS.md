# V8 exact constructions and their limits

These are supplied mathematical constructions, not acquired research methods.
The physical premise remains a finite, time-independent Lindblad generator.
All approximation certificates are conditional on that generator describing the
intended system; none repairs model error by numerical exactness.

## 1. Overlapping anticommuting coefficient covers

Let R = Σ_p r_p P_p be a Hermitian Pauli expansion with real coefficients. For
pairwise anticommuting sets G, choose x_Gp with exact reconstruction

    r_p = Σ_{G containing p} x_Gp.

Set A_G = Σ_p x_Gp P_p. Since P_p²=I and distinct anticommuting pairs cancel,

    A_G² = (Σ_p x_Gp²) I,
    ||R||∞ ≤ Σ_G ||A_G||∞ = Σ_G sqrt(Σ_p x_Gp²).

This proof does not require positive splits, disjoint sets, an optimizer, or
optimality. An accepted rational upper u_G≥0 with u_G²≥Σ_p x_Gp² suffices.
The independent checker re-adds every coefficient, checks every pair, checks
each squared upper bound, and sums the advertised bound exactly.

The proposer may use floating arithmetic. It rescales coefficients by their
largest magnitude, proposes allocations, rationalizes them, and places the
exact discrepancy for each coefficient into one of its incidences. The
checker never trusts the floating objective or a convergence claim. Arithmetic
and parser bounds apply to the submitted rationals; a larger problem may be
refused even when a mathematical cover exists.

### A strict separation from every disjoint partition

Take the five two-qubit Paulis IX, IY, XI, YX, ZY, each with coefficient 1.
Their anticommutation edges are

    IX—IY, IX—ZY, IY—YX, XI—YX, XI—ZY.

This graph is a five-cycle, so there is no clique of size three. The best
partition contains two edges and one singleton, with value 2√2+1. Assigning
half of each coefficient to each incident edge instead gives five edges of
norm √(1/2), for total 5/√2. The improvement is strict because

    2√2+1 − 5/√2 = 1−1/√2 > 0.

It is approximately 7.65% relative to the optimal disjoint partition bound.
This is a separation of certificate forms, not a demonstrated reduction in
complete computational cost. On the complete three-vertex clique X,Y,Z,
the full partition bound √3 beats the three half-edge cover 3/√2; the baseline
must receive that full clique. A previous agent test missed this comparison
and was discarded.

### A coordinate minimization identity

Hold every split except those for one coefficient c fixed. Let

    a_j = sqrt(Σ_{q != p} x_jq²),     Σ_j x_jp = c.

The triangle inequality in R² gives

    Σ_j sqrt(a_j²+x_jp²) ≥ sqrt((Σ_j a_j)²+c²).

If a=Σ_j a_j>0, equality is attained by x_jp=c a_j/a. This proves a global
minimizer for that single coordinate block. If a=0, any same-sign split has
objective |c|. This elementary fact justifies the floating coordinate proposal;
it does not prove global convergence of cyclic minimization, or monotonicity
after finite rationalization. Only the resulting feasible certificate counts.
The implementation also retains conventional partition candidates and returns
the smaller verified bound. Group generation and every sweep are charged.

## 2. Exact residual certification under the extended norm interface

For O(t)=exp(t L†)O0 and a Hermitian polynomial approximation P(t), define
R=P'−L†P. Positive unital Heisenberg propagation contracts the operator norm
on Hermitian matrices, so variation of constants yields

    ||O(t)−P(t)||∞ ≤ ||O0−P(0)||∞ + ∫_0^t ||R(s)||∞ ds.

For piecewise polynomials, add the norm of each actual interface jump. A
continuous numerical trajectory has zero jumps even though its earlier
residual error remains in the accumulated integral. Previous truncation error
must not be counted again as an interface jump, and must not be discarded.

The extended V8 checker recomputes the same exact residuals as V7, but accepts
coefficient-cover witnesses for their norm. Partition and overlap arms use the
same checker. All coefficients, initial differences, horizon equality, and
claimed sums are independently checked. The cover checker permits up to 512
groups so the V7 512-term singleton partition remains admissible; the search
family remains capped at 128 groups and 2048 incidences. The first timing run
exposed the interface mismatch at 128 checker groups and was preserved as an
invalidated comparison before the corrected replay.

## 3. Clearing denominators in the exact Taylor recurrence

Let d clear all rational Hamiltonian and damping coefficients, so D=d L† has
integer Pauli columns. Let q0 clear the initial observable coefficients and
b0=q0 O0 be integral. Define

    b_(k+1)=D b_k,      c_k=b_k/(q0 d^k k!).

Induction gives c_(k+1)=L†c_k/(k+1). Therefore P_m(t)=Σ_(k=0)^m c_k t^k is
exactly the ordinary Taylor polynomial, not an approximation to the previous
implementation. Its residual is −L†c_m t^m. Python integer accumulation avoids
repeated rational reduction inside multiply-adds, while exact Fraction exports,
group construction, bit bounds, cache bounds, and independent residual checking
remain charged. The full exported coefficients were compared across every
successful development case, not just final predictions.

This is conventional denominator clearing. Its positive timing result improves
the supplied baseline; it is not evidence that Spectra acquired m1.

## 4. A nonlinear parity of the operator graph

Every Hamiltonian term in the current development families contains an even
number of Y factors, hence H is a real matrix. A Hermitian Pauli word is real
or purely imaginary according to the parity of its Y count. Multiplication by
i in i[H,P] swaps these sectors. Thus the Hamiltonian commutator graph is
bipartite with charge

    c(P)=#Y(P) mod 2 = Σ_i x_i z_i mod 2

in binary symplectic coordinates. This charge is quadratic, explaining why a
search restricted to linear symplectic charges missed it. Local depolarization
is diagonal and preserves the charge, so the damped generator as a whole is
not an off-diagonal bipartite operator. The existing sparse action already
omits forbidden edges; the symmetry alone does not demonstrate fewer operations.

## Method-acquisition gate

These results establish mathematical validity or conventional speed, not the
requested causal chain. A future m1 must be an acquired executable construction
rule that transfers, and its availability must reduce the fully charged cost
of acquiring a useful m2. Final S0/S1/S2 increments must pass on untouched
systems. Cost per candidate and number of candidates are separate outcomes.
Amortization may span a declared downstream horizon, but acquisition and
validation costs cannot disappear from that ledger.
