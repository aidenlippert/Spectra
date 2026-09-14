# Matched M=12, N=6 Reynolds quartic certificate

This is a bounded size-transfer experiment for the matched half-filled model
with 12 spin orbitals, 6 particles, hopping `t=1/5`, and degree-4 SOS terms.
The Reynolds construction averages over flavor permutations and global
flavor swap, then expands the averaged positive Gram matrices before exact
certificate replay.

## Structural reduction and timings

The unreduced degree-4 construction has 1,145 charge blocks. The symmetry
calculation retains 30 charge orbits and 216 reduced PSD blocks. The largest
reduced block has dimension 14 and the total number of symmetric PSD scalar
variables is 2,599. The invariant coefficient equations reduce from 3,640
rows to 39, with multiplier rank 19.

Assembly took 20.503008 s. Clarabel solve time was 0.489196 s and returned
`optimal`. Export and exact interval replay are included in the receipt and
completed successfully.

## Certified interval

The exact replay gives the lower endpoint

`2687792842876988729/500000000000000000 = 5.375585685753977`.

The independently generated Rayleigh upper endpoint is

`26877938273878065/4999999911102389 = 5.375587750350995`.

Thus the certified width is `2.0645970172874137e-6`. The raw SOS residual

`1032422511271/500000000000000000 = 2.064845022542e-6`

is charged by the exact checker; no floating-point objective is treated as a
certificate.

For the matched model, the independent sector reference identifies the
ground sector as zero doubly occupied pairs and brackets its energy by
5.375587750350633 and 5.375587750351133 (width `5e-13`). The certificate's
upper endpoint is slightly below the reference upper endpoint, so the two
intervals do not nest; they overlap, with the reference value lying inside
the certificate interval. The exact checker also reconstructs the
Hamiltonian and confirms it equals the matched collective Hamiltonian.

## Interpretation

The M=12 run shows that the Reynolds reduction remains small at this next
size: 2,599 reduced PSD variables versus 2,082 for the M=10 run under the
same reduced formulation. (The much larger M=10 unreduced count is not an
apples-to-apples comparison.) It also produces a valid, independently replayable
energy interval. This is evidence for transfer of the symmetry reduction to
one larger matched instance; it does not establish asymptotic scaling,
accuracy for generic Hamiltonians, or a universal physical-marginal
representation.
