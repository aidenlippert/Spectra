# A changed proposal geometry with the same exact energy certificates

The Hamiltonian, ten spatial patterns, collective tail and entire allowed
SOS family are unchanged. Every energy certificate is accepted by the same
[spin-completion checker](/Users/aidenlippert/Documents/Spectra/research/spin_completion_20260913/core.py).
The existing exact CAR, number-ideal, residual and tail argument in the
[preceding proof](/Users/aidenlippert/Documents/Spectra/research/spin_enrichment_20260913/PROOF.md)
applies. No new energy-validity assumption is introduced by the pricing rule.

## The operator-size metric

For each fixed symmetry frame Phi, define G_y,ij = y({Phi_i^dagger,Phi_j})
and T_ij = tau({Phi_i^dagger,Phi_j}). Here tau is the normalized uniform
trace on the full N-electron sector. Its diagonal occupation moments follow
the existing exact combinatorial seed function; no sector enumeration is
needed. For B = sum_i v_i Phi_i, the quotient

    v^T G_y v / (v^T T v)

measures the proposed violation relative to tau(B^dagger B + B B^dagger).
The denominator is a uniform-sector mean squared operator size. It is not
a ground-state expectation, a new physical measurement or a wavefunction.

The generalized eigenproblem G_y v = lambda T v proposes its extreme
directions. This quotient is invariant under invertible changes of frame
coordinates when both Grams transform together. Ordinary eigenvalues using
the Euclidean norm of frame coefficients do not have that invariance.
Neither choice guarantees the best energy gain; exact candidate solves test
that question.

## Stable numerical realization

Some coefficient-frame trace matrices have condition estimates around
3.6e10. Direct generalized eigensolves failed the coordinate-rescaling test;
diagonal equilibration alone did not fix it. The final implementation works
before the ill-conditioned coefficient contraction.

Let C contain the exact operators' coefficients converted to floats in a
canonical monomial basis, and compute the full thin QR decomposition C=QR.
For the monomial moment matrices M_y and M_tau, form

    G_Q = Q^T M_y Q,    T_Q = Q^T M_tau Q.

Solve G_Q z = lambda T_Q z, with invertible diagonal equilibration, then
recover v from R v = z. No column is dropped, no ridge is added, and a rank-
deficient frame or non-positive trace metric is refused. This is a numerical
proposal calculation only. The energy solver's separate conditioning is
unchanged; the QR transform does not densify its full-frame SDP maps.

Before rounding, v is rescaled to Euclidean length one, preserving its
direction. The rounded operator column C v is used to recompute both the
ordinary and trace-normalized moments without subtracting nearly cancelling
entries of the ill-conditioned coefficient Gram. Independence is checked in
the same operator coefficient space as before. The recorded negative
threshold is now in normalized units, so the tested proposal rule includes
both its metric and that normalized eligibility threshold.

The final root and direction factors are still rationalized and expanded
exactly. A wrong proposal, numerical eigensolve, rounding effect or solve
error cannot silently improve the bound: the reconstructed coefficient-L1
residual penalty is charged by the unchanged verifier.

## What a comparison can establish

The source certificate, dual hashes, Hamiltonian, full frame, baseline and
360-second budget must match the historical control. All completed valid
new candidates count, including those not selected. The best verified
interval at a resource ceiling is the minimum among actually completed
certificates within that ceiling, including the inherited source. The
comparison does not interpolate missing runs or replace a stronger earlier
certificate with a weaker later result.

The new certificates, inherited source and old control best are freshly
replayed using only the standard library. Other old checkpoints reuse their
preceding exact replays with matching certificate hashes. One such historical
comparison does not establish general superiority, statistical reproducibility
or scaling. The 1.6 mHa target and small-representation question remain
separate from a relative improvement in either curve.

A conditional full-family dual diagnostic uses the preceding exact checker
and repair protocol. Only an accepted full-family witness plus the separately
replayed precise physical proof can establish a new true energy-error floor.
A failed or weak repair does not rule out this operator family.
