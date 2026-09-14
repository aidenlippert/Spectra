# V4: finite Clifford interaction family for learned experiment design

Use (n) qubits and candidate Pauli interactions (P) of weight at most two.
A primitive gate is

\[
U_P=e^{-i\pi P/4}.
\]

An episode contains either one gate or an ordered pair (U_{P_2}U_{P_1}),
drawn from a declared finite candidate family. The learner receives experiment
records but no hidden support or circuit-family labels.

## Exact likelihood

Prepare the physically separable Pauli-biased state

\[
\rho_{s,O}=\frac{I+sO}{2^n},\qquad s\in\{\pm1\},
\]

where (O) is a Pauli string. Measure Pauli (Q) after the circuit (C). In
the Heisenberg picture compute (C^\dagger Q C=\sigma R), where (R) is a
Pauli string and (sigma\in\{\pm1,\pm i\}) is the Pauli multiplication phase.
For a valid Hermitian observable, the surviving phase is real and

\[
\mu(s,O,C,Q)=\operatorname{Tr}(\rho_{s,O}C^\dagger Q C)
 =\begin{cases}s\sigma,&R=O,\\0,&R\ne O.\end{cases}
\]

An ideal (+1) outcome has probability

\[
p_+=(1+\mu)/2\in\{0,1/2,1\}.
\]

With a declared visibility/readout channel ν represented by a rational
visibility (v\in[0,1]), use

\[
p_+=(1+v\mu)/2.
\]

Thus (N) shots have an exact binomial likelihood

\[
\Pr(k|e)=\binom Nk p_+^k(1-p_+)^{N-k},
\]

with rational (p_+) whenever (v) is rational. The likelihood must include
channel noise and calibration errors; an unexplained empirical Bernoulli rate is
not a physical source model.

## Resource ledger and policy task

An experiment record includes ((s,O,C,Q,N)), outcomes, preparation cost
(w(O)), measurement cost (w(Q)), gate count (|C|), and channel-noise cost.
The policy chooses the next experiment by expected information gain per total
cost, using the finite candidate likelihoods. A dynamic-programming policy over
the same finite family provides a conditional optimum for validation.

Task variants must change more than coefficient scale: use an unseen candidate
interaction family, an ordered-pair circuit family, and a new coupling topology.
Also vary public channel noise and its cost. A policy that learned a real rule
should change measurement and preparation choices when noise makes a formerly
useful measurement less informative; a fixed schedule is a negative control.
Evaluate on withheld circuits and interactions, with frozen, retrieval-only,
from-scratch Bayesian, random, and exhaustive finite-family baselines. Match
shots, gate calls, preparation weight, and computation.

## Operational assumptions

The source assumes separable preparation of ((I+sO)/2^n), exact Pauli control,
known Clifford phase conventions, and a calibrated rational noise channel.
Preparation is not free: weight and sign-selection costs are charged. The
candidate family is finite and public to the evaluator; otherwise finite-grid
posterior guarantees do not apply. Identical likelihoods under all allowed
experiments are an equivalence class and must trigger abstention. Agreement on
the finite family cannot establish a law outside it.

The construction is a clean bridge from a learned experimental strategy to a
physical likelihood because every prediction reduces to Pauli conjugation and a
finite binomial model. It still establishes only conditional transfer and
resource accounting, not open-ended mechanism discovery.
