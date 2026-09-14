# Direct-discovery audit contract

This audit applies to fixed sparse-atom LP and locality-block SDP proposals.
A smaller final certificate is not sufficient: the result must show that the
retained atoms were found from the Hamiltonian and declared metadata, and that
the omitted part has a sound, independently replayed error bound.

## Required evidence

1. **Source independence.** Discovery may read only the Hamiltonian, sector
   labels, declared structural metadata, and its own generated state. It must
   not read an FCI vector, an existing certificate's factors, lower endpoint,
   upper witness, or a hidden full dictionary. A provenance test should inject
   false source metadata and verify identical discovery output or explicit
   refusal.
2. **Exact atom validity.** Every retained factor must replay as a rational
   positive square. For LP atoms, the coefficient map must be the exact CAR map
   with the correct dagger convention. Quantized coefficients must be replayed
   after rounding; floating solver residuals are diagnostics only.
3. **Omitted residual coverage.** The residual must be formed from the full
   Hamiltonian identity, including number-ideal terms and every higher-degree
   word. Any high-degree word not handled by a locality or block theorem must
   be charged by a declared exact norm bound. A coefficient residual over only
   degree-4 rows is insufficient when generated products have higher degree.
4. **LP/SDP signs and sectors.** Check the primal identity and objective signs
   on a one-atom hand case. The multiplier ideal must vanish on the declared
   fixed- N sector, and every Gram block must have a definite charge/parity
   sector. A numerical solver status is not an exact feasibility certificate.
5. **Locality-block soundness.** A block norm or interface penalty must be
   independently proved for the stated sector. Local minima from floating
   eigensolvers need outward-rounded lower bounds before they can enter a
   claimed global interval. Overlap ownership must not double count or leave
   an uncovered interaction.
6. **Total cost.** Report candidate generation, rejected candidates, LP/SDP
   solves, exact replay, upper witness generation, maximum block dimension,
   and all determinant/state/operator enumeration. A compact output that first
   scans an exponential dictionary does not pass discovery scaling.
7. **Held-out transfer.** Freeze the selection rule before applying it to new
   geometry, coupling, basis, or active space. Same-geometry representation
   changes are useful robustness checks but do not count as held-out chemistry.

## Current expected verdict

The method passes only if the exact replay validates the complete emitted
certificate and the error bound is no larger than the reported interval. It
may still be useful if it fails scaling or transfer, but the report must call
that a bounded heuristic or a post-processing result. No direct-discovery
experiment currently establishes universal efficient discovery or generic
FeMoco capability.

## First implementation audit

`direct_sparse_discovery.py` currently fails the direct-discovery contract.
Its equality matrix includes `+I` and `-I` slack columns, then computes
`eq @ sol.x - rhs`; that quantity is identically zero by construction and is
not the physical coefficient residual. The physical residual must be formed
from the identity columns only, with positive and negative slack values
charged separately. The run therefore cannot report `residual_l1=0` as an
SOS result. It also mixes annihilation words of charge -1 and -2 in paired
atoms, while the production certificate schema requires each Gram block to
have one charge. Finally, it has no certificate export or exact replay, and
invoking the file directly fails because it does not establish the repository
module path; the module invocation works. Its `atom_support` records only the
first word and loses the second word and sign, so provenance is incomplete.

`direct_local_grams.py` builds blocks from Hamiltonian supports and does not
read source factors or upper witnesses for construction. A smoke run on the
M4 certificate produced a production-checker-accepted export, but its receipt
is still a numerical SDP proposal (`proposed_b` and floating solver state).
The exact checker acceptance is the authoritative fact; the solver residual
and reported lower are not independently sufficient. The locality block rule
also needs a declared omitted-block/error comparison before it can claim
scaling: the current run prices all generated local blocks and does not bound
the contribution of ungenerated nonlocal or high-degree words.

## Recheck of the rewritten sparse method

The rewritten `direct_sparse_discovery.py` fixes the earlier defects. It builds
the complete residual word set, including generated degree-six terms where
present, exports rational factors, and calls the production `verify()` checker.
Independent diagonal and off-diagonal four-mode smoke cases produced valid
rational lower certificates; the off-diagonal toy reached the exact lower
value `-1`. The `m4_a128` artifact has no independent upper witness, so it is
a valid lower construction rather than a complete energy interval.

The method still scans 2,828 word pairs for the H4 molecular runs and 14,718
for H6 before selecting a bounded atom set. This cost is recorded, and
`omitted_family_optimality_proved` is false. That flag does not weaken lower
bound soundness because the full residual is charged; an omitted-family oracle
is unnecessary for validity. It does mean efficient discovery scaling remains
unproved.

## Exact dual obstruction audit

`python -S research/certificate_scaling/direct_dual_replay.py` passes. The
replay checks `y_1=1`, `||y||_infty<=1`, all 443 body-2-or-lower number-ideal
constraints, 18,496 singleton and signed equal-magnitude pair-square
constraints, and computes the Hamiltonian pairing from exact source words.
It returns `U=-4.459007306545172 Ha`, source lower `L=-3.330388761351986 Ha`,
and `L-U=1.128618545193186 Ha`.

The sign is correct for this restricted cone: for every represented identity
`H=bI+ideal+positive_squares+R`, dual feasibility gives
`b-||R||_1 <= y(H)=U`. The independent source lower therefore quantifies the
gap between the known ground-energy lower evidence and the best lower endpoint
available to this dictionary/residual objective. Cross-charge pairs are
included. The replay also finds a negative determinant for one unrestricted
two-word PSD block, showing that the tested singleton/signed-pair cone is
strictly narrower than arbitrary pair Gram PSD.

This is an exact obstruction for the declared H4 quadratic dictionary with a
body-2 ideal and full coefficient-l1 residual. It is not an obstruction to
all SOS certificates, arbitrary PSD pair blocks, higher-degree atoms, or the
physical ground energy itself. The Lambda interval replays and this dual
obstruction answer different questions and should not be merged into a generic
scalability claim.
