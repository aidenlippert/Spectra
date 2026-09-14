# Adaptive target-only lower bound for a hidden Pauli prefix

## Model and calibration conditions

Fix r = j - 1 hidden action axes A = (a1,...,ar), uniform on a menu of size 2^r. Let b be an independent uniform target bit. A target-only shot chooses an action q from the same menu, possibly adaptively from all previous actions, outcomes, and private randomness. Conditional on a prefix match q = A, the shot has a known positive signal of visibility v: its outcome distribution depends on b. Conditional on q != A, the outcome is a fair coin independent of b. The hidden A and b are sampled once per target instance.

The result assumes prior calibration has not selected or biased A: A remains uniform conditional on the learner's public prior information. Target shots reset to the same source state, so a nonmatching shot does not alter later hidden axes or latent state. A stateful apparatus is a different model; if a failed query changes A, the state, or later outcome laws, the theorem must be redone for that transition kernel. Calibration observations that leak A, or a nonuniform posterior after calibration, replace 2^-r by the largest posterior mass of a query and may remove the bound.

## Coupling theorem

Let H_t be the event that one of the first t queries matches A. Let N be a null process with the same adaptive policy, but every target outcome is an independent fair coin. Couple the real and null processes by using the same policy randomness and the same fair outcome whenever the real process has not yet matched A.

**Theorem.** Under the reset model and uniform independent A, for every adaptive t-shot policy,

Pr(H_t) <= min(1, t / 2^r).

Moreover, if P0 and P1 are the transcript laws conditional on b=0 and b=1, then their total variation distance is at most t / 2^r. Consequently the Bayes error for recovering b obeys

BayesError >= 1/2 - t / (2 * 2^r).

The error bound is deliberately loose: it claims only a sign/target-bit advantage of at most t / 2^(r+1), and does not claim exact zero error at any finite t.

**Proof.** Construct one coupling for P0 and P1 using the same hidden A, policy randomness, and independent fair null outcomes. Run both conditional processes on the shared null path until its first query matching A. Before that hit, every real outcome is fair for either b, so the two transcripts coincide and their adaptive queries remain the same. The first hit in either real process occurs exactly when the shared null path first queries A; after that point the coupling may diverge. A is independent of the entire null path, including its adaptively generated queries, so a union bound over the at most t null queries gives Pr(hit) <= t 2^-r, capped at one. Hence TV(P0,P1) <= Pr(hit). The binary-testing identity for equal priors gives Bayes error = (1-TV(P0,P1))/2, proving the bound. This avoids conditioning on “no hit,” under which A would no longer be uniform. □

The coupling also shows why adaptivity does not defeat the bound: adaptivity changes q_s as a function of the fair null history, but that entire null path is independent of A. The argument does not claim that the real transcript remains null after a hit.

## Tightness and small-menu check

For t <= 2^r, an open-loop policy that queries t distinct prefixes has hit probability exactly t / 2^r. Thus the probability bound is tight under the stated menu. Exhaustive enumeration of distinct-query policies gives:

| r | menu size | maximum hit probability for t = 1,2,3,4 |
|---|---:|---|
| 1 | 2 | 1/2, 1 |
| 2 | 4 | 1/4, 1/2, 3/4, 1 |
| 3 | 8 | 1/8, 1/4, 3/8, 1/2 |

The adaptive policy cannot improve these values, because no pre-hit transcript distinguishes the hidden prefixes. Repeating a query is weaker than enumerating new prefixes.

## Signal visibility and finite shots

The lower bound above concerns the event of ever seeing a signal. If a matching shot has visibility v < 1, even a known match has single-shot Bayes error (1-v)/2 for a symmetric binary sign channel. Conditional on H_t, later inference can be optimal and still cannot beat the event-level bound. A simple mixture argument therefore gives the displayed error lower bound independently of v; keeping v explicit only strengthens practical bounds.

If a policy receives k independent shots per queried prefix, count them as k target shots: the union bound becomes Pr(H) <= k t / 2^r when t denotes distinct query attempts, or <= T / 2^r for total shots T. A matching query may then be repeated to estimate its sign, but repetition does not increase the probability of discovering the hidden prefix.

## Boundaries and falsification tests

The theorem does not cover calibration that changes the prior, a stateful source, correlated hidden axes, actions that reveal partial prefix information while nonmatching, or a measurement whose outcome is biased even on mismatch. Each is an identifiable model change. The correct diagnostic is to estimate the null conditional law for every action before first match and test whether it is independent of A and b.

The result is a finite prior-and-menu lower bound, not an unrestricted statement about quantum experiments. It proves that target-only adaptation cannot replace calibration when the hidden prefix is uniform and all nonmatching observations are exactly uninformative. A conventional Bayesian policy with the same prior can attain the same optimal menu strategy; no learner superiority follows from this theorem.
