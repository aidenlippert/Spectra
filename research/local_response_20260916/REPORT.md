# Local construction, exact family limit, and collective correction

**The requested field-level many-body breakthrough remains unsolved.** This
pass built and exactly replayed a compact lower-bound component on a fully
interacting eight-site Hubbard ladder. It also proved why two tested levels
of approximation cannot deliver the desired accuracy. It did not establish
controlled recursive cost, a new tractable class, or new thermodynamic physics.

This is the open 2-by-4 Hubbard ladder at U/t=8, with four up and four down
electrons. It is not the molecular hydrogen H8 fixture. Its energies are in t
units, not hartrees; all previous molecular accuracy milestones are preserved.

## What was constructed

Three overlapping four-site blocks include every rung and leg hopping, every
onsite interaction, and all local charge sectors. Arbitrary number-preserving
one-rung boundary operators cancel between adjacent blocks. A direct local
semidefinite optimization proposes those operators. An independent rational
checker verifies the actual matrices and every positivity assertion.

An additional collective correction uses local spectral profiles and the
incompatibility of overlapping low-energy subspaces. The accepting proof uses
at most 36-by-36 local PSD matrices, 16-by-16 reduced local matrices, and a
scalar projector-overlap inequality. Neither lower construction uses a global
determinant basis or a global ground-state teacher.

| Accepted quantity | Value in t units | Meaning |
|---|---:|---|
| Zero-boundary local lower | -4.4925060891 | Matched local baseline |
| Optimized boundary lower | -4.3751736540 | Exact physical lower component |
| Boundary-family ceiling | -4.3751532484263764 | Ceiling on **every** lower certificate in the declared family |
| Additional collective correction | +0.2678136631 | Exact nonlocal compatibility penalty |
| Corrected lower | -4.1073599909 | Stronger exact physical lower component |
| Separately certified variational upper | -3.025922805690684… | Uses 4,900 enumerated global configurations |
| Complete interval width | 1.0814371852093165 | Misses the predeclared 0.001t target |

The word “nonlocal” in the correction describes the globally valid operator
inequality. Its proof does not construct a full global matrix. It contracts
local data and uses exact projector algebra.

## What is now proved insufficient

The boundary-family optimum lies between the accepted lower and its exact
locally consistent marginal witness, a bracket only 0.0000204055736236t wide.
This is not an optimizer-stall diagnosis. That entire declared local family
has a proved ceiling.

The collective correction exceeds the old family ceiling, giving

\[
E_0-D_{\mathrm{old}}\ge 0.2677932575263764t.
\]

This proves a physical shortfall of the old family without assuming that a
numerical reference eigenvalue is exact. The correction genuinely adds
information absent from the old family.

The rank-one spectral compression itself is also insufficient. A locally
specified trial state for the *coarsened profile Hamiltonian* gives an exact
upper ceiling of about 0.568t on any universal correction that follows from
only the three retained local ground projectors and gap inequalities. The
full local Hamiltonians contain additional excitation information. This is
an information limit for that compression, not an upper bound on the physical
Hubbard ground energy and not an impossibility theorem for richer methods.
The exact value and matched-upper implications are in `RESULT.json`.

A numerical screen retaining ranks 1,2,4,6,8,12,16 in the same grouped-projector
rule did not outperform rank one. All those attempts are recorded. They do
not rule out other rank choices, multilevel spectral profiles, or different
collective response constructions. At rank 16 the overlap contraction already
contained 1,048,576 entries, illustrating why rank growth must be charged.

## What the proof actually does

With exact positive local blocks A_i, certify

\[
A_i\succeq g_i(I-P_i).
\]

The outer projectors P_0,P_2 commute and their product is rank one globally.
Set g=min(g_0,g_2), h=g_1. A local contraction gives
p=||(P_0P_2)P_1||^2. Then

\[
\sum_i A_i\succeq\gamma I
\]

whenever the exact rational checks establish

\[
0\le\gamma\le\min(g,h),\qquad(g-\gamma)(h-\gamma)\ge ghp.
\]

The full derivation, the dual proof, and the limits of the construction are in
[DERIVATION.md](DERIVATION.md). These use established projection and
quantum-marginal techniques. Local SDP lower bounds already have primary
literature precedents, including [Eisert](https://doi.org/10.1103/cz6k-y46r)
and [Fawzi, Fawzi, Scalet](https://arxiv.org/abs/2311.18706). We do not claim
that this pass discovered a new general principle or established priority.

## Verification, dependencies, and cost

Eight focused tests cover both physical hopping legs and spin signs, global
model reconstruction, all-charge local coverage, boundary-message symmetry,
dual consistency, the projector contraction against an independent small
matrix, the profile-model witness against a sparse physical state, and
invalid-profile rejection. Four altered certificates are rejected by fresh
replay. A regression enumerates the 4,900-state sector specifically to check
that the local representation reconstructs the original model exactly; that
is validation work, not the accepting local checker.

The final fresh replay runs with Python's standard library only. It rechecks
the primal, dual, collective correction, profile information limit, and a
separately enumerated physical upper. Its accepted interval still fails the
accuracy target. Solver “Solved” flags are never acceptance criteria.

The local boundary constructor took 2.99 seconds including exact primal/dual
repair, peaking near 259 MB; its optimization used 55 scalar variables and
75 PSD blocks. The rank-one correction constructor took 1.47 seconds. These
are separately measured component processes, not a fresh complete-solver
speedup claim. The complete accounting includes three timeouts: two stalled
CVXPY import probes and an unrelated import chain in an initial checker test.
The final implementation calls Clarabel directly and uses a standalone exact
local PSD routine. No paid compute, installation, commit, or push was used.

The ledger is `results/local_response_20260916/ACCOUNTING.json`. Reasoning,
code authoring/review, literature work, early unmetered agent checks, and report
assembly are excluded from subprocess totals. Agent-proposed errors in the
norm bound, model hopping, spin validation, and source metadata were corrected
before acceptance. Exploratory agent notes are not substitutes for the
consolidated derivation and actual replay.

## What remains missing

The earlier response theorem's retained positivity and relative residual
metric have not been closed at useful accuracy. This pass supplies a compact
alternative global correction and an exact obstruction to its first two
simplifications. It does not prove that a richer retained excitation spectrum
and collective response can be constructed, contracted, and verified at
controlled cost through repeated elimination.

There is no new many-body solution, tight enumeration-free full interval,
scaling theorem, or competitive physical prediction from this pass. The new
result is an exact local-to-global improvement with precisely measured limits.

## Replay

Run from the Spectra root, with a new output directory:

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -B -S -m research.local_response_20260916.replay \
  --certificate results/local_response_20260916/boundary/certificate.json \
  --correction results/local_response_20260916/projector/correction.json \
  --upper results/local_response_20260916/reference/upper.json \
  --out results/local_response_20260916/another_fresh_replay
```

Omit `--upper` to replay only the local lower component, its family obstruction,
and the collective correction without global enumeration. Exact endpoint
strings in `results/local_response_20260916/final_replay/receipt.json` are
authoritative; decimals in this report are display values.
