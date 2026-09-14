# Frontier routes (57–64)

This report covers quantum-assisted energy estimation, symbolic/formal
discovery, learned proposal systems, and exact fallbacks for certified
molecular ground-energy intervals. It deliberately separates a genuine
certificate from a heuristic approximation.

## Concrete result: overlap barrier

`python research/all_angles_20260913/frontier/overlap_barrier.py` gives:

```
p=1/4: E[shots]=4, shots_90=9
p=1/16: E[shots]=16, shots_90=36
p=1/100: E[shots]=100, shots_90=230
```

For the two-level Hamiltonian `H=diag(0, Delta)` and trial state
`|phi>=sqrt(p)|0>+sqrt(1-p)|1>`, ideal phase estimation reports the ground
energy with probability exactly `p`. Independent repetitions therefore have
geometric mean `1/p`, and at least `ceil(log(1-alpha)/log(1-p))` shots for
success probability `alpha`. This is an exact lemma, checked with rational
arithmetic. Amplitude amplification can improve query scaling to
`O(1/sqrt(p))`, but requires coherent state preparation, reflections, and a
known/isolated target window. Thus a classical GPU simulation cannot claim this
quantum speedup.

## Route status and falsifiable tests

| Route | Candidate for Spectra | What would count as success | Main blocker |
|---|---|---|---|
| 57 Phase estimation | Prepare a symmetry-adapted trial, simulate `e^{-itH}`, phase-estimate energy, then convert an upper estimate into a lower certificate using a separately certified gap/residual bound | Interval width and total state-preparation/simulation cost beat current exact/SDP baseline on a fixture | Fault-tolerant Hamiltonian simulation; overlap `p`; phase-estimation output alone is not a lower bound |
| 58 Quantum Krylov | Measure `S_ij=<phi_i|phi_j>` and `H_ij` for `{A_i|phi>}` and solve generalized Ritz problem; certify omitted-space residual by block inequality | Same interval target with fewer operator terms than direct certificate, with interval arithmetic for measured entries | Noisy matrix elements, ill-conditioned `S`, and omitted-space bound |
| 59 Filtering/amplification | QSVT/polynomial filter around a guessed ground window, then phase estimation; use classical certificate to validate window | Reduced effective overlap cost on fixtures with known spectral gap | Requires block encoding and gap; filter leakage must be bounded |
| 60 Analog/quantum simulation | Use analog simulator only to propose orbitals/ansatz; replay proposal classically and certify algebraically | Better candidate state while certificate remains classical and reproducible | Device calibration and no rigorous energy guarantee without replay |
| 61 Symbolic identity discovery | Search commutator/anticommutator identities and SOS regroupings over a small operator grammar; exact rational replay in SymPy | New positive decomposition with fewer monomials/Gram rank than baseline | Expression search grows combinatorially; numerical fits are not proofs |
| 62 Formal lemma search | Emit Lean/Isabelle proof obligations for positivity, CAR reductions, Schur complements, and interval arithmetic; accept only kernel-checked proofs | A machine-checked certificate generated from a discovered identity | Formalization cost and missing library lemmas |
| 63 Learned proposal policy | Train/rank proposals from Hamiltonian graph, orbital integrals, and failed certificates; every proposal is exact-replay checked | Higher verified-certificate yield per CPU/GPU hour on held-out molecules | Distribution shift and reward hacking; learned score has no evidentiary value alone |
| 64 Exact fallback | Sparse/full diagonalization, interval Sturm/inertia, branch-and-bound over occupation sectors; use as small-fixture oracle | Definitive interval or obstruction and regression oracle | Exponential scaling; useful for calibration and falsification only |

## Proposed hybrid experiment

For each fixture, generate a classical trial state (selected-CI or tensor
network), use Krylov vectors to estimate a low-dimensional Ritz value, and ask
the certificate engine for a complement lower bound. If the complement bound
is `QHQ >= mu Q`, the Schur bound is

`E0 >= theta - ||QH phi||^2/(mu-theta)` for `mu > theta`, while
`E0 <= theta` is variational. The exact two-level check above shows why the
quantum route must report the trial overlap and gap assumptions explicitly.

The decisive benchmark should record (a) certificate width, (b) discovery
time, (c) verification time, (d) number of terms/Gram variables, (e) trial
overlap or residual, and (f) assumptions needed. A GPU run can accelerate
linear algebra or candidate search; it cannot substitute for a quantum phase
estimation oracle or make an uncertified approximate energy rigorous.

## Sources (primary)

- Quantum Krylov ground/excited energy algorithms: https://arxiv.org/abs/2109.06868
- Exact quantum Lanczos construction: https://quantum-journal.org/papers/q-2023-05-23-1018/
- Shorter phase estimation for ground energies (Ding–Lin): https://arxiv.org/abs/2211.11973
- Ground-state preparation and overlap scaling: https://arxiv.org/abs/1607.05625
- Formal verified mathematics workflow (Lean): https://doi.org/10.1145/3371078

These references motivate routes; none establishes a Spectra certificate or a
new complexity theorem. No novelty claim is made.

## Route-64 exact active-space oracle

`research/all_angles_20260913/frontier/exact_h4_oracle.py` reconstructs the
assigned active-space-ladder H4 fixture at the absolute source path
`results/certificate_scaling/active_space_ladder/h4/fixture.json`. Its SHA256 is
`8e64505deb4cc06a87f75e669b0f3073d15ebca72daa6d086e06eeea0b0f5120`.
The 70-dimensional fixed-N sector receives an exact rational LDL PSD lower
certificate and a rationalized numerical eigenvector upper replay:

* lower `-3.66699997`, upper `-3.6669999559641995`;
* certified width `1.40358006675e-8` Hartree;
* actual exact work is cubic in dimension with large rational intermediates,
  while building the sector is exponential in modes/particles.

The earlier `-4.4759` result is preserved as a legacy STO-3G stress fixture,
with SHA256 `1157834626d04c6a71e90ddc785b66440ca5472a735918c66b0efa578a21bdf2`;
it is not the active-space result. Tests reject a downward corruption, reject
the same-dimension legacy fixture by content hash, and include a zero-pivot
coupled indefinite matrix regression. The exact checker validates square,
Hermitian matrices and binds modes, particles, and fixture content.
