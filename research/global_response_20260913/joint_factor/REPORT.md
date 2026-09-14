# Existing joint-factor audit

The prior cubic certificate was checked against the H6 fixture used by the
global-response work.  The canonical Hamiltonian hashes match exactly:

`14e7229dd5001f634e352dd92346f2f4a67792003789cb4abe9e36a79d2a88`.

The independent spectral replay gives width `0.000054994392835... Ha`, with
12, 66, and 220 dimensional exterior-power components for bodies 1, 2, and 3.
The lower replay enumerates no many-body determinants.  This is already a
complete compact-factor lower proof for the frozen H6 Hamiltonian and is much
stronger than the 1.6 mHa target.

The qualification is more precise: the earlier cubic discovery itself uses
only the Hamiltonian and sector data, and does not enumerate determinants.
What it does not provide is a cheap *small-factor* discovery rule: its full
cubic dictionary and numerical solve are expensive, and the resulting factor
description remains large. The next meaningful global-response experiment is
therefore to construct analogous factors directly from the response-reduced
operator and compare exact interval and cost against this baseline.

As a bounded direct-family control, the existing Hamiltonian-ranked two-word
program was run from frozen H6/H8 Hamiltonian inputs only. H6 budgets 64 and
256 produced exact lower candidates -42.6448 and -41.3853 Ha, respectively;
H8 budget 64 produced -90.6440 Ha. Their residual L1 penalties dominate, so
none is a useful energy bound. They are valid negative controls: merely
selecting a small deterministic family does not reproduce the full cubic
factor. The runs constructed no Fock basis and used no stored factor or
amplitude inputs; timings were 0.52 s, 1.85 s, and 1.65 s.
