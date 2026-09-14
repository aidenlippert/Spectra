# V4 finite Bayesian experiment-policy certificate

V4 separates an exact teacher from a policy learned from data. The teacher is
allowed to solve a bounded finite dynamic program; the learned tree is evaluated
independently on the same task and is not allowed to reuse the teacher's action
labels at inference time.

## Fixed finite task

The input contains at most eight hidden models `m`, eight actions `a`, two
outcomes `y`, a rational prior `p(m)`, rational likelihoods `Q(y|m,a)`, a
rational per-action cost, a shot cost, and a horizon `H<=3`. Every table row is
included in the certificate and hashed. Rows must be nonnegative and normalize
exactly: `sum_y Q(y|m,a)=1`; priors normalize exactly as well.

The posterior is represented by an integer weight vector `w` obtained by
clearing the common denominator of the prior and likelihood table, then
canonicalizing by the gcd. This is a finite exact state representation. A
history is never identified only by a floating-point feature vector.

## Exact teacher recursion

Let `R(w)` be Bayes sign/classification risk at the current posterior, or a
declared terminal loss over models. With remaining depth `h` and remaining
resource budget `b`, define

`T(w,h,b)=R(w)` at `h=0`, and otherwise

`T(w,h,b)=min_a [ cost(a) + sum_y P(y|w,a) T(w_{a,y},h-1,b-cost(a)) ]`.

An action is admissible only when its cost is within `b`. If action costs are
not charged as an objective, they must still be reported and bounded; otherwise
the teacher can select arbitrarily expensive experiments. For each action the
certificate stores every outcome probability, canonical successor posterior,
successor value, total expected cost, and value. The verifier checks all
admissible actions, not merely the selected one, and checks that the selected
value equals the exact minimum. A state-count limit of 5,000 causes explicit
abstention rather than silent truncation.

The teacher certificate proves optimality only for this finite table and budget.
Memoization is a computational optimization, not an extra source of knowledge:
the verifier reconstructs every reachable state and rejects dangling nodes.

## Learned policy and posterior features

The learned policy is a finite decision tree whose internal node stores a
declared posterior feature vector and an action label. A safe feature is an
exact, task-defined function of canonical weights, such as normalized model
masses, entropy represented by an indexed lookup cell, or likelihood-ratio
interval bins. Arbitrary neural embeddings are not certificates because equal
features need not imply equal future risk.

To claim generalization, define finite feature cells by exact predicates, for
example `w_i/w_total <= q` or membership in a listed rational interval. The
certificate must prove that all test histories assigned to a cell use the same
action and report cell coverage. The learned tree is evaluated by exhaustive
Bayesian rollout over the supplied model table, not by replaying only the
training histories. Its risk and expected resource cost are independently
summed as rational numbers.

The comparison reports teacher risk/cost, learned-tree risk/cost, regret, and
the fraction of reachable posterior mass covered by declared cells. An
uncovered cell is an abstention or failure, never silently assigned the teacher
action.

## Fairness and anti-oracle controls

Training and teacher planning use the same model, prior, likelihood, action,
noise, and cost tables. The teacher may use exhaustive dynamic programming only
within the declared state and horizon limits. The learned policy may use data
generated from the same task but may not receive teacher action labels unless
the experiment explicitly studies imitation learning. If teacher labels are
used, label-generation compute and sample cost are reported separately.

Posterior memoization is allowed for both systems as a runtime optimization.
Caching a precomputed optimal action table is not allowed in the learned policy
unless that table is itself the learned artifact under evaluation. Hidden
planning oracles, uncharged extra shots, changing priors, and adaptive changes
to the likelihood table invalidate the comparison.

## Required certificate checks

The verifier rejects malformed rational values, non-normalized priors or
likelihoods, unsupported actions, budget overruns, duplicate or unreachable
nodes, noncanonical posteriors, missing outcome branches, and learned nodes
whose action is not defined for the declared feature cell. It independently
recomputes the teacher Bellman minimum and learned-policy Bayes rollout.

The certificate should include a JSON summary with table hashes, horizon,
budget, state count, teacher and learned values, expected experiment cost,
regret, coverage, and an explicit `abstained` flag. A positive result means the
learned policy transfers within the finite declared task. It does not establish
physical-model correctness or unrestricted scientific discovery.
