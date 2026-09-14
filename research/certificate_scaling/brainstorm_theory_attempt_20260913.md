# Brainstorm: trial guided certificates and coupled energy bounds

This note responds to the proposed “complete landscape” as a research map. I
agree with its engineering direction—adaptive local structure, explicit
residuals, and an independently replayable one sided certificate—but not with
any implication that one route is known to scale generically. The strongest
claim currently supported by the repository is a conditional program: search
for a compact restricted certificate, then certify the omitted directions by a
separate operator bound. Existing notes already cover Galerkin/Arnoldi,
Duhamel references, wedge/Gershgorin bounds, commutator dictionaries, cluster
tails, and locality obstructions. The ideas below add combinations around a
trial state and a coupled energy difference.

## One proved modest lemma (Schur/Kato style)

Let a self-adjoint finite matrix be written relative to a normalized trial
vector `u` and its orthogonal complement as

```text
H = [ theta   r* ]
    [ r       C  ],       theta=<u,Hu>,  r=QHu.
```

Assume the independently certified complement inequality `C >= mu I`, with
`mu > theta`. Then

```text
lambda_min(H) >= theta - ||r||^2/(mu-theta).
```

Proof: set `L = theta - ||r||^2/(mu-theta)`. The lower-right block of `H-LI`
is at least `(mu-L)I`, and its Schur complement is bounded below by
`theta-L-||r||^2/(mu-L)`. Writing `d=mu-theta` and `x=||r||^2`, this is
`x/d - x/(d+x/d) >= 0`. Hence `H-LI >= 0`. This is an actual certificate
provided `theta`, `||r||^2`, and `mu` have outward-rounded rational/algebraic
bounds. It is the same basic mechanism behind Kato–Temple residual estimates,
but stated in the form directly useful for a certificate checker.

The attractive failure mode is also clear: a tensor-network parent Hamiltonian
gap can supply `mu` only when the complement is genuinely the orthogonal
complement of a certified ground space (or a certified low-energy subspace).
An approximate MPS/DMRG vector does not supply that fact. Injectivity can imply
a gap for a parent Hamiltonian, but the physical `H` need not equal that parent
Hamiltonian; perturbation transfer requires a norm bound and a gap margin.
Likewise, a gap estimate for `QHQ` cannot be inferred from a small trial
residual. The cross block `r` is known and enters the Schur correction, but a
small residual alone does not certify that `u` is the ground state: for
`H=diag(0,1)` and `u=e2`, the residual is zero while the complement block is
`0 < theta=1`. This is the proof boundary, not an implementation detail.

A clean sufficient parent transfer condition is: `P_parent u=0`, `u` is the
unique parent ground vector, `Delta>0`, `alpha>=0`, `Q P_parent Q >= Delta Q`,
and `H=alpha P_parent+V` with `Q V Q >= v_Q Q`. Then the needed complement
bound is explicitly `mu=alpha Delta+v_Q`. Establishing these inequalities for
the physical model, including the correct `Q`, is the hard part.

## Routes worth trying

1. **Trial-state-guided SOS with Schur completion.** Obtain a low-rank trial
   state (MPS, selected CI, or Krylov vector), form its exact local energy and
   residual, and ask SOS only for a lower bound on the complement block. The
   checker combines the SOS bound `mu` with the lemma above. The hypothesis is
   that the complement bound is much easier than a global SOS. It is unproved;
   test first on small Hubbard chains where the full spectrum is available.

2. **Joint interval optimization.** Optimize a trial upper state and a dual
   lower certificate jointly, with a penalty or explicit constraint on their
   residual coupling. Use `E_trial - E_lower` as the target while retaining a
   separately valid lower endpoint at every iterate. This could make operator
   selection focus on directions that affect the current interval width. It is
   only safe if each iterate remains a valid one-sided certificate; a small
   observed gap is not itself a proof. A distinct, more speculative variant is
   a coupled ground-energy *difference* certificate for two Hamiltonians, where
   shared trial structure and a certified bound on `H_1-H_2` might cancel
   extensive errors; no such cancellation is assumed here.

3. **Parent-Hamiltonian bridge plus perturbation accounting.** Construct a
   frustration-free local parent `P` for a trial tensor network and write
   `H=P+V`. Certify `P>=0`, a finite-size gap for `P`, and `||V||` or a sharper
   blockwise bound. Weyl/Davis–Kahan style perturbation then controls the
   low-energy sector, after which the Schur lemma handles the trial residual.
   This is theoretically clean but likely loses its advantage for long-range
   Coulomb terms or a poor tensor ansatz.

4. **Residual-weighted local SOS / sparse elimination.** Rank omitted operator
   blocks by their contribution to `QHu` and by a Schur sensitivity estimate
   `||r_X||^2/(mu-theta)`, then solve local SOS certificates in batches. Couple
   this with the existing wedge residual and cluster-tail bounds so that blocks
   with no local mass are eliminated analytically. The new conjecture is a
   safe, computable stopping rule based on the sum of block sensitivities; it
   must be checked against full SDP improvements on toy systems.

5. **Hybrid analytic/numerical resolvent.** Use a certified reference
   propagator or Krylov compression for the trial sector and SOS only for the
   resolvent/complement. This combines the repository’s Duhamel and Galerkin
   ideas with the Schur formula. The main risk is nonnormality and coefficient
   growth; no generic contraction or polynomial scaling follows.

Other theoretical avenues remain: quantum belief propagation and quasi-adiabatic
continuation, Lieb–Robinson cluster bounds, stoquastic Perron bounds,
representation-theoretic block diagonalization, interval verified Lanczos,
noncommutative Positivstellensatz hierarchies, and randomized trace/residual
sketches followed by deterministic replay. These are candidate tools, not
evidence that the target problem is solved.

## Falsifiable next experiment

For `H` on a 4–8 site chain, choose an MPS or Krylov trial `u`. Compute an exact
rational enclosure for `theta` and `x=||QHu||²`. Independently produce a
complement lower certificate `mu` (start with the exact restricted matrix, then
try local SOS/parent-Hamiltonian bounds). Compare the Schur endpoint with the
true smallest eigenvalue and with the existing global certificate at equal
coefficient/checker budget. The route fails its intended purpose if `mu<=theta`,
if the parent-to-physical perturbation consumes the gap, or if the endpoint is
consistently weaker than the baseline after accounting for trial-state and
complement-certificate work.

## Agent prompt

> On a small fixed Hubbard/XXZ chain, investigate the block certificate
> `H=[theta,r*;r,C]`. Produce an independently replayable rational lower bound
> on `C`, combine it with `theta` and `||r||²` through the Schur inequality, and
> compare against exact diagonalization and the existing certificate. Attempt
> both a tensor-network parent-Hamiltonian complement and a local SOS bound.
> Record every assumption, failed case, coefficient/checker cost, and whether
> the physical-vs-parent perturbation leaves a positive gap. Do not claim a
> generic gap theorem or external “Navier–Stokes agent” access.

## Primary sources

- [Kuroda, Kato–Temple estimates (RIMS, 2007)](https://www.kurims.kyoto-u.ac.jp/~prims/pdf/43-2/43-2-24.pdf) — residual/eigenvalue inequality.
- [Saad, *Numerical Methods for Large Eigenvalue Problems*, 2nd ed.](https://www-users.cse.umn.edu/~saad/eig_book_2ndEd.pdf) — approximate eigenvectors and Kato–Temple exposition.
- [Fannes, Nachtergaele & Werner, Finitely correlated states (Commun. Math. Phys., 1992)](https://doi.org/10.1007/BF02099414) — parent-Hamiltonian/MPS gap foundations.
- [Nachtergaele, The spectral gap for some spin chains with discrete symmetry breaking (Commun. Math. Phys., 1996)](https://doi.org/10.1007/BF02101844) — martingale gap method.
- [Rozmán, Molnár & Schuch, Lower Bounds on Spectral Gaps of Parent Hamiltonians via Tensor Networks (2026)](https://arxiv.org/abs/2607.19078) — recent computable parent-gap refinements.
- [Putinar, Positive polynomials on compact semi-algebraic sets (Indiana Univ. Math. J., 1993)](https://doi.org/10.1512/iumj.1993.42.42017) — SOS/Positivstellensatz basis.
- [Lasserre, Global optimization with polynomials and the problem of moments (SIAM J. Optim., 2001)](https://doi.org/10.1137/S1052623400366802) — SDP moment/SOS hierarchy.
- [Bravyi, Hastings & Verstraete, Lieb–Robinson bounds and the generation of correlations (Phys. Rev. Lett., 2006)](https://doi.org/10.1103/PhysRevLett.97.050401) — locality/propagation premise for cluster-style bridges.

Novelty has not been checked against the complete literature. None of these
references establishes the proposed scalable chemistry certificate.
