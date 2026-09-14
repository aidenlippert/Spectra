# Schur trial-state route: bounded executed attempt

## Scope and method

This is a small exact rational experiment of the trial-state Schur route. The
code uses a basis-vector trial (a controlled analogue of a selected-CI/HF
determinant), computes `theta=<u,Hu>`, the exact residual square
`||QHu||^2`, and a rational complement lower bound by Gershgorin. If
`mu>theta`, it emits `theta-||r||^2/(mu-theta)`; otherwise it refuses. A NumPy
eigensolve is used only as a post-freeze validation oracle. No ground-state
eigenvector is used to construct a trial.

## Results

The parent-like 3x3 rational example had `theta=0`, residual square `1/100`,
and complement bound `mu=9/5`, producing the certified endpoint `-1/180`.
The oracle minimum was approximately `-0.0050208119`, so the endpoint is below
it as required. The naive norm endpoint `theta-||r||` was `-0.1`, substantially
weaker here.

The deliberate excited-eigenvector case `diag(0,1)` with `u=e2` had zero
residual but `mu=0<=theta=1`, so the implementation refused. This confirms
that residual zero alone cannot identify the ground state. A third 3x3 case had
`mu=1/2<=theta=1` and also refused; its naive norm endpoint was `0.2928932`,
above the true oracle minimum `0.5` only because the trial was not a valid
ground-state lower-bound construction. The Schur gate correctly avoids this
unsupported inference.

## Original rational molecular H4 fixture

The saved straight H4 STO-3G fixture was reconstructed directly from its 70
determinants (8 spin orbitals, 4 particles) using the repository CAR action.
The trial was determinant bitmask `15`, a Hamiltonian-derived RHF candidate
also present in the independent upper-state artifact; no FCI eigenvector was
used. The exact Rayleigh value was
`-440714970861/125000000000` (about `-3.5257198`) and the exact residual square
was `23777795598926675102523/250000000000000000000000` (about `0.0951112`).
Full 69x69 complement Gershgorin gave
`mu=-4155084279457/1000000000000` (about `-4.1550843`), so `mu<=theta` and the
Schur certificate correctly refused. The cheaper diagonal-only value was
`-3405168600593/1000000000000`; this is diagnostic only, not a valid complement
bound, because it ignores off-diagonal couplings. The exact oracle
minimum was approximately `-3.666999956`; the saved independent upper energy
is linked in `h4_run.json`. Full complement work was 4,761 matrix entries,
small here but explicitly exponential in general determinant-space dimension.

## Parent/local route assessment

The full-complement Gershgorin bound is an oracle-scale baseline: it is easy to
replay but requires a complement matrix whose size grows with the full Hilbert
space. A parent route would need a certified `P_parent u=0`, unique parent
ground vector, `Delta>0`, `alpha>=0`, `Q P_parent Q>=Delta Q`, and
`H=alpha P_parent+V` with `QVQ>=v_QQ`; then `mu=alpha Delta+v_Q`. No physical
Hubbard/molecular parent transfer was established in this run. Naive global
perturbation (`||V||`) was not competitive conceptually: it gives a useful
bound only after a certified parent gap and often consumes the entire margin.

## Route status ledger

The route numbering refers to the requested all-angles map; only statuses are
claimed here, not prior-art or novelty.

| Route | Status in this bounded run |
|---|---|
| 21 | Executed as exact Schur endpoint on a rational trial; accepted on one case. |
| 22 | Executed as residual-square input; no independent scalable complement bound. |
| 23 | Refusal gate tested on an excited eigenvector; residual alone rejected. |
| 24 | Full-complement Gershgorin used as oracle baseline; scaling unresolved. |
| 33 | Parent-Hamiltonian transfer condition derived; no physical instance certified. |
| 34 | Local/parent complement route remains unimplemented hypothesis. |
| 35 | Naive norm perturbation comparison recorded; weaker or margin-consuming. |
| 36 | Molecular H4 target not attempted in this five-minute bounded preflight. |

## Reproduction and limits

Run the command in `receipt.json`; the focused unittest has two passing tests.
The experiment is exact in its rational arithmetic for `theta`, residual, and
Gershgorin/Schur values. The floating eigensolve is validation only. It does
not establish a chemistry result, scalable discovery, or a tensor-network gap
for the physical Hamiltonian. The next meaningful experiment is an H4 or small
Hubbard determinant-space trial with a genuinely local parent inequality and
an independently checked `QVQ` lower bound.
