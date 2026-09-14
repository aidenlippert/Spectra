> Proposal status: unexecuted, not an established V5 result. The initial full-matrix parameter-count advantage was rejected; equally informed baselines can exploit the admitted operator family. Later qualifications are essential.

# V5 black-box joint-design test

## Exact finite construction

Each episode samples a hidden state (s\in\{0,1\}^d) uniformly. An action chooses a subset (S\subseteq[d]) and receives a noisy bit

\[
y=\bigoplus_{i\in S}s_{\pi(i)}\oplus e,
\qquad e\sim\mathrm{Bernoulli}(\eta),
\]

where (pi) is an unknown permutation and the target label is a known Boolean function of (s) supplied after the episode. The learner sees only ((S,y,	ext{target})) and does not receive a likelihood table or the hidden state. A calibration action may also choose a known reference state, allowing later readout-flip calibration.

The finite hypothesis class is a declared permutation family (Pi) and parity supports of size at most (q). All action costs, episode counts, (eta<1/2), and confidence (delta) are fixed before testing. Equivalence classes that induce identical action distributions are intentionally indistinguishable.

## Two successive discoveries

Discovery A tests whether univariate action statistics predict the target. If every single-coordinate action is at chance but a two-coordinate action has nonzero correlation, the learner must introduce a joint feature. It then chooses pair actions by information gain over the surviving permutation/support hypotheses. The held-out test samples a new parity support and scores target prediction and the number of episodes required to identify it.

Discovery B occurs only after A: use the newly identified joint support to select reference actions whose hidden parity is known, and estimate a readout-flip parameter or coordinate permutation ambiguity. The planner then chooses the cheapest discriminating action for a second held-out support. This is a change in experiment selection and source calibration, not merely another regression record.

The positive criterion is sequential: after A, held-out experiment cost falls relative to a univariate-only planner; after B, cost falls again on a new support with the same hidden source. Costs include every exploratory action, failed hypothesis test, reference shot, and target observation.

## Strong baselines and accounting

Use the same records for four systems: (i) univariate-only learner; (ii) passive all-pairs learner; (iii) generic Boolean learner with the declared permutation/parity class; and (iv) stateful active planner that can select actions but is forbidden from storing a discovered reusable relation. The generic learner gets the same hypothesis class and records; any gain must come from the learned representation and action policy, not privileged metadata.

Shared hidden-state episodes require joint-source accounting. Reusing an episode across candidate actions is invalid unless the physical protocol permits simultaneous readout; otherwise each action consumes a fresh episode. Correlated outcomes within an episode require a concentration argument over independent episodes, not independent-shot Hoeffding. A readout calibration must be reported as a separate source and cost.

## What this does and does not establish

This is closely related to V2's finite GF(2) response-law learner: both identify a parity law from noisy binary observations. The additional test is active experiment design under shared latent episodes, with a second discovery that changes which calibration action is cheapest. It does not establish open-vocabulary physical discovery, and a generic Boolean learner may match it. A result counts only if the representation improves held-out action choice after amortizing discovery and beats the matched active generic baseline.

Cheap falsifiers are: univariate correlations already reveal the target; passive all-pairs search has lower total cost; the second calibration does not alter the optimal action; or the gain disappears when every reused episode is charged correctly. These outcomes are useful negative results because they show V4 regression improvement was record fitting rather than reusable joint structure.
