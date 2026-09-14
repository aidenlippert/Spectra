# Structural certificate research — first goal campaign, 2026-09-12

The active goal is to find a useful, efficiently checkable structural condition
that controls both certificate discovery and verification as systems grow.
The chemistry scaling breakthrough remains open. We now have stronger molecular
controls, an exact counterexample to an unjustified locality inference, and
working scalable positive controls for two restricted noncommuting families.

Continuation: `HIDDEN_BASIS_RESULTS.md` records Hamiltonian-only recovery of
hidden rational orbital structure through 16 modes, and an exact quantitative
obstruction on the molecular ladder. The obstruction covers arbitrary complex
orbital rotations and one-body number-ideal multipliers, within the stated
quartic coefficient residual norm. It does not obstruct general SOS factors.

## Corrections to the proposed landscape

- Chordal PSD completion, sparse PSD decomposition, and a sparse fermionic SOS
  relaxation are different statements. A molecular orbital graph does not
  automatically give a valid low-width Gram decomposition. Sparse
  noncommutative hierarchies already exist; convergence alone is not a bound on
  the degree needed for fixed accuracy. See [Klep–Magron–Povh](https://arxiv.org/abs/1909.00569)
  and [Wang–Magron](https://arxiv.org/abs/2010.06956).
- A gap plus suitable locality can imply ground-state clustering. It does not
  imply that a relaxed pseudomoment dual obeys the same decay. This distinction
  matters even for a connected path with a uniformly gapped vacuum. See
  [Hastings–Koma](https://arxiv.org/abs/math-ph/0507008) and the exact counterexample below.
- Exact symmetry can safely remove cross-sector Gram variables; it does not
  prove bounded sector multiplicities. This campaign uses verified Abelian
  alpha-charge and GF(2) parity symmetries, not a newly implemented SU(2) quotient.
- `python -S` demonstrates standard-library independence. Polynomial replay
  requires bounded factor degree and accounting for certificate size and
  rational bit lengths. It is not a physical-system-size guarantee by itself.
- Four increasing active spaces are a useful falsification ladder, not an
  asymptotic theorem. The target is fixed TOTAL interval width, here 0.0016 Ha,
  not an error allowance independently granted to each site or interaction.

## Molecular ladder: real execution, remaining failure

Generated separate linear H4/H6/H8/H10 chains at 1.4 Å spacing, STO-3G,
neutral total electron count, in RHF canonical and occupied/virtual-block
localized bases. These are distinct from the earlier square/rectangle H4
fixtures. Each Hamiltonian is rationalized and its coefficient difference
from the floating integral construction is recorded. This is not a certified
continuum or integral-evaluation error.

Twenty quadratic/adaptive/control jobs completed on the two existing A10
hosts; four additional cubic-ladder jobs produced three certificates and one
300-second H10 timeout with no exported certificate. All 23 remote certificates
were independently replayed with exact sparse upper witnesses. A subsequent
fixed-factor number-ideal repair of H6 also replayed. All 95 downloaded result
file hashes matched. Remote CPU execution is reported; no GPU speedup is claimed.

Best canonical ladder intervals in this bounded campaign:

| Chain | Best returned method | Certified width (Ha) | Target |
|---|---|---:|---|
| H4 | Full cubic, exact symmetry | 0.0000486761 | Pass |
| H6 | Full cubic + ideal repair | 0.0068454618 | Fail |
| H8 | Full quadratic, exact symmetry | 0.0269978844 | Fail |
| H10 | Full quadratic, exact symmetry | 0.0768817055 | Fail |

The H8/H10 upper witnesses retain 1,000 FCI coefficients. Their truncation
error contributes to interval width. FCI generation is explicitly a validation
cost, not a scalable upper-state discovery result. Upper calculations use
integrals reconstructed from the frozen rational Hamiltonian; an initial
regenerated-orbital attempt was misaligned and is superseded by the
`active_space_ladder_references_aligned*` outputs. No upper witness or old
factor certificate enters lower-certificate discovery.

The old H4 square and rectangle controls also pass with new full cubic
symmetry-reduced certificates: widths 0.0000448304 and 0.0000079900 Ha,
respectively, at 97,122 and 78,868 bytes including upper witnesses. These are
finite benchmark results, not a new complexity or prior-art claim.

The uniform five-round adaptive-small-block rule deteriorates badly on larger
chains. This rejects that fixed search budget as a working scaling strategy;
it does not prove that every sparse block representation fails.

Exact symmetry reduced canonical H10 quadratic Gram entries from 233,000 to
44,300. Full-SDP H6 symmetry/no-symmetry runs took approximately 1.22/7.87
seconds and gave matching numerical bounds up to solver error. Single-run
comparisons are not repeat timing statistics.

## Graph route: a concrete obstruction for this representation

All eight tested orbital graphs are complete, in both canonical and block-
localized coordinates. Their exact treewidths are 3,5,7,9 along the ladder.
Even after coefficient pruning charged a full 0.0016 Ha L1 norm budget, every
edge remains. The diagnostic also computes the coefficient-L1 cost of deleting
ALL terms inducing each edge. Every minimum edge-deletion cost exceeds the
budget; for localized H10 it is about 1.435 Ha. Thus arbitrary coefficient
removal under this budget cannot remove even one edge of this particular
orbital graph. This does not exclude a different operator graph, norm bound,
quotient, interaction factorization, or orbital transform.

Receipt: `results/certificate_scaling/structural_graph_diagnostics.json`.

## Exact counterexample to gap-only dual pruning

For H=sum n_i - (1/4)sum nearest hopping on a path, the full-Fock vacuum has a
uniform gap at least 1/2 and zero disjoint connected correlations. A restricted
LP retaining individual monomial squares, with no number ideal, admits an
optimal dual whose omitted endpoint moment block is [[0,1],[1,0]]. Its normalized
negative Rayleigh quotient is -1 at distances 3,7,11,15. The exact CAR replayer
checks all retained squares, normalization, residual-dual limits, and matching
primal/dual objectives. It claims a pricing violation, not an objective gain
from adding that endpoint atom. It does not refute a theorem for a stronger
number-ideal relaxation.

Use `locality_dual_counterexample.py` and `locality_dual_exact.json`; earlier
agent-only concatenation/perturbation diagnostics are superseded.

## Additional route 1: positive amplitudes plus bounded local-energy width

For a real stoquastic H and strictly positive psi, define local energy
L(x)=(H psi)(x)/psi(x). The ground energy lies between min L and max L.
When log psi is a sum of bounded-scope features and the Hamiltonian moves have
bounded scope, L can be a sum of local functions. A supplied bounded-width
tree decomposition permits dynamic-programming extrema, without global
configuration enumeration.

This has a direct positivity bridge: each negative off-diagonal transition
with rate t gives a PSD 2x2 block t[[psi_y/psi_x,-1],[-1,psi_x/psi_y]]. Local
ratios permit local conditional blocks. DP message cancellation certifies the
remaining diagonal lower bound. Rational factor weights permit exact replay.
The log-parameter hypograph is convex, but a polynomial discovery theorem also
needs a certified separation oracle, bounded parameters/bit cost and a useful
approximation family. Our scalar numerical optimizer is not that general theorem.

Executed control: periodic spin chains with positive Jastrow amplitudes
q^(-domain walls), q inferred by optimization from H only. These are deliberately
constructed solvable, noncommuting models. At 8,16,32,64,128,256 spins every
returned interval is exactly [0,0]. The 256-spin certificate is 10,900 bytes;
discovery took 0.535 s and replay 0.027 s. Exact DP work is linear in length
at fixed rational bit complexity; no global state vector is constructed.
Independent small dense-matrix, brute-force DP, perturbation and rejection
controls pass. This is a restricted spin-chain backend, not the molecular CAR
solver. See `stoquastic_chain_certificate.py` and `stoquastic_chain_audit.md`.

## Additional route 2: compatible local kernels in an interacting fermion chain

For positive rational-square r, consider the open-chain Hamiltonian

H = sum_i [r n_i + r^(-1)n_(i+1) - (r+r^(-1))n_i n_(i+1)
           - a_i^dagger a_(i+1) - a_(i+1)^dagger a_i].

Each term is exactly Q_i^dagger Q_i, where
Q_i=sqrt(r)a_i(1-n_(i+1)) - r^(-1/2)a_(i+1)(1-n_i).
The factors have degree three; the Hamiltonian is quartic and adjacent energy
terms do not commute. In every fixed-N sector, amplitudes proportional to
r^(sum_i i n_i) are annihilated by every local term. Their norm is computed by
a particle-count DP. No many-body state list is needed.

The compiler infers r from H's boundary coefficient, checks the ENTIRE
Hamiltonian pattern, constructs factors and proves a matched upper bound. It
rejects a coefficient mutation or mismatched upper parameter. This is an
H-only recognizer of a known solvable exclusion/ferromagnet-like family,
not a discovered solution of general molecular Hamiltonians. Under the usual
open-chain fermion/spin correspondence this is a ferromagnetic XXZ chain
with telescoping boundary fields; the corresponding solvable kink-state
structure is established prior art ([Koma–Nachtergaele](https://arxiv.org/abs/cond-mat/9709208)).

For M=8,16,32,64,128,256 orbitals at half filling, all intervals are exactly
[0,0]. At M=256,N=128: 255 factors, 1,020 factor coefficients, 81,788 bytes,
0.068 s discovery, approximately 1.34 s replay including the exact norm DP.
The norm DP has O(MN) updates and polynomially growing integer bit lengths;
only factor storage is claimed linear. Independent small full-sector action
checks cover r=1,4,1/4 and all N=0,...,4.

See `fermionic_ratio_chain.py`, its audit, and
`results/certificate_scaling/fermionic_ratio_chain/summary.json`.
This positive control is not completion of the chemistry-oriented goal.

## Other routes and next target

- **Weak-coupling linked clusters with a certified tail.** If the aggregated
  omitted tail is at most N A exp(-nu k), taking k=O(log(N/epsilon)) can maintain
  fixed total error with polynomial cluster work under bounded branching.
  Ground-energy algorithms in this regime already exist
  ([Bravyi–DiVincenzo–Loss](https://arxiv.org/abs/0707.1894)). The missing step
  here is an efficiently constructed positive normal form and certified
  uniform tail, not just a perturbative energy series. Cluster counts must be
  included in nu; local coefficient decay alone is insufficient.
- **Recover a hidden interaction basis.** Analyze the kernel of the one-body
  commutator map A -> [V,Q(A)] for the quartic interaction V. A maximal
  commuting self-adjoint algebra with simple joint spectrum can identify a
  density-density orbital basis. Then test graph width and local-ratio
  conditions in that recovered basis. Dimension alone is insufficient;
  commutativity and rotation error must be checked. This is a next experiment,
  not an implemented result in this campaign.
- **Low-rank interactions, approximate Markov recovery, and impurity structure.**
  These are useful representations or restricted-family hypotheses; each still
  needs a one-sided certificate and an aggregate error bound. Low tensor rank
  or an accurate MPS alone does not imply a compact lower certificate.

The next useful advance is extending a proved local positive structure to
less manufactured Hamiltonians, including recovery from unknown orbital
coordinates and certified perturbations at fixed interaction strength. Simply
increasing sparse-atom counts or presenting the solvable controls as generic
chemistry would not satisfy the active goal.
