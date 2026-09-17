# Constructive response investigation

Date: 2026-09-16. This pass tests a proposed research direction. It does not
assume that a world-level breakthrough can be delivered on demand.

## Fixed objective

Construct and rigorously bound the eliminated response of an interacting
fermion Hamiltonian without a hidden full-state teacher, and preserve a useful
representation and error bound through repeated elimination. Established Schur,
Lanczos, perturbative renormalization, and variational residual identities are
baselines, not new discoveries.

## Independent investigations

- Exact residual bounds, positive preconditioners, and recursive composition.
- A primary-source novelty audit, including matrix resolvent quadrature and
  provable impurity algorithms.
- Hubbard local coercivity with charge transfer and overlap accounting.
- Counterexamples to unjustified gap/locality/compression implications.
- Existing molecular implementation and accepting-checker audit.
- Direct experiments on fully interacting Hubbard ladders, comparing bare and
  locally dressed retained sectors of the **same** Hamiltonian.

The parent owns this protocol, numerical Hubbard diagnostics, integration,
accounting and final report. Agents have separately assigned files. Existing
certificates and completed campaigns are not modified.

## Computational controls

Use the repulsive Hubbard model at U/t=8 on declared open 2-by-L ladders.
Vertical dimer hopping is t=1; horizontal hopping lambda is continued from
0 through 1/4, 1/2 to 1. Lambda=1 is a fully interacting endpoint. Half filling
uses the complete N_up=N_down=L sector; no fragment charge is fixed.

The initial locally dressed projector is a product of exactly orthogonal
rational two-site rotations selected solely from the isolated local bond.
Its angle is not fitted to a global ground state. Bare and dressed calculations
use exactly equivalent global Hamiltonians, apart from explicitly distinguished
floating proposal arithmetic. Full diagonalizations are charged diagnostics,
not claimed non-enumerative discovery or accepted molecular certificates.

The disconnected dimer family is an analytic structural control only. It
cannot satisfy the interacting-system objective.

Run only one sizeable local numerical process at a time. Record failed and
timed-out attempts, elapsed time and peak child memory. Use bounded subprocesses
and existing libraries. No paid compute, package installation, Git writes or
GitHub operations are part of this pass.

## Acceptance and reporting

Accept algebraic results only with checked hypotheses and exact independent
small examples or counterexamples. Distinguish numerical diagnostics from
exact rational replay and conditionally rounded arithmetic. Count determinant
enumeration wherever it occurs. A compact program or exported file is not a
complexity result.

For a claimed breakthrough require a new, precisely scoped result beyond the
identified prior art, a complete constructive/verification argument and actual
evidence for its significant physical or computational scope. Otherwise report
the proved component, precise obstruction, and unsolved obligation honestly.
