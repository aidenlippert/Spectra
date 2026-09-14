# Overlapping cluster lower certificates

`overlap_cluster_lower.py` assigns each original Hamiltonian CAR word to the
first cluster containing its support (weight one; hence weights sum exactly),
builds the complete local Fock matrix, and certifies its minimum with exact
rational Gershgorin bounds over all local particle sectors. Nonfitting words
are charged by exact coefficient L1. No global Fock enumeration is used.

| variant | assigned/total | omitted L1 | global lower |
|---|---:|---:|---:|
| H4 disjoint 4 | 28/184 | 15.8473 | -21.5932 |
| H4 overlapping 6 | 108/184 | 6.8311 | -12.4797 |
| H6 disjoint 4 | 66/918 | 52.2399 | -63.0256 |
| H6 overlapping 6 | 166/918 | 37.3515 | -50.0460 |

The exact rational values, block dimensions, and minima are in
`results/wave2_20260913/geometry/overlap_cluster_lower_receipt.json`.
The table rounds only the compact summary; exact values are authoritative.

The overlapping windows reduce the charged residual on both fixtures, but the
corrected H6 disjoint-4 comparison has three genuine width-4 blocks (the
previous run accidentally used an 8-mode block). Gershgorin is still loose.
The next structural improvement is
an exact rational LDL/interval eigenvalue certificate for each local block;
the assignment and residual accounting can remain unchanged. This experiment
does not assume telescoping cancellation or infer treewidth from block size.
