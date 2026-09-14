# A finite adaptive-experiment cost lower bound

Use a two-stage binary sensor with a reusable nuisance mode. The hidden mode is
`m in {0,1}` and the fresh target is `s in {-1,+1}`. The sensor action `x`
returns

`y = s` when `m=0`, and `y=-s` when `m=1`.

All variables and outcomes are finite and deterministic; a rational noisy
variant replaces the equality by a binary symmetric channel with crossover
`q=1/4`. A calibration action `c` reveals `m` exactly at cost one. A target
read costs one. The learner is not given the mode label as a source annotation:
it must obtain it through `c` or infer it from labelled calibration outcomes.

## Two successive tasks

In stage 1, the system receives a calibration target whose sign is known. It
may spend one read on `c`, after which it stores the reusable representation
`mode -> sensor polarity`. In stage 2, a fresh independent sign is presented
under the same mode, with no new mode label.

With the retained mode representation, stage 2 needs one target read and has
zero error. A frozen/scratch policy that discards stage-1 information has three
options:

* one target read: its posterior on `s` remains uniform over the unknown mode,
  so error is `1/2`;
* two target reads: repeated deterministic reads do not identify the mode, so
  error remains `1/2`;
* one calibration read plus one target read: cost two, zero error.

Thus at zero-error threshold the exact experimental-cost gap is one read per
fresh target task: retained knowledge costs one target read, while scratch
requires calibration plus target read. The result is a nuisance-identification
gain, not a new parity-computation example.

For the noisy `q=1/4` channel, after mode calibration one target read has Bayes
error `1/4`. Without calibration, the mixture of the two opposite channels is
exactly uniform, so every number of target reads has posterior error `1/2`.
The same strict gap holds for any target-error threshold below `1/2`.

## Exact certificate and fairness boundary

Represent the task as four finite hypotheses `(m,s)` with rational prior
`1/4`, two target outcomes, and a calibration outcome table. A finite Bayes DP
enumerates both actions and proves the minimum expected cost for a declared
error threshold. The retained policy's state includes the mode-polarity
mapping; scratch resets it before stage 2. Both policies use identical target
priors, likelihoods, action menu, and measurement costs.

The “same-data conventional learner” can tie the retained policy's predictions
on stage-2 outputs if it is allowed to retain the calibration mapping. That is
the correct baseline: the theorem demonstrates cumulative reuse relative to a
frozen or retrained-from-scratch policy, not an algorithmic superiority theorem
against a stateful conventional implementation with the same memory and data.
Acquisition and storage costs must be reported separately. A representation
that stores a raw calibration example rather than the polarity rule has not
shown the claimed reusable operation until it transfers to a changed target
prior or a new sensor instance.

The obstruction is equally clear: if the mode changes independently before
each stage, the retained mapping has no value and the strict savings vanish.
Therefore the compositional claim is conditional on a declared invariant
nuisance parameter, and the mode-change experiment is a direct falsifier.
