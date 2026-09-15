# Accepted development results

These are conditional mechanism tests. They are not fresh integral-to-certificate
timings. Complete cold comparisons and held-out results are reported separately.

## A direct H10 proof meets the original accuracy target

The local-orbital construction with overlapping windows of widths 2, 3 and 4,
plus one coherent collective dictionary containing pair-supported cubics and
the declared three/four-orbital cubics, certified the **original H10 Hamiltonian
to 1.5565392773062308 mHa**. The independent checker reconstructed the rotation,
both energy endpoints, spin allowances and all residual coefficients.

The constructor generated **702,332 Gram entries**, with largest block 422,
instead of the preserved canonical reference's **1,287,700 entries**, largest
block 505. This is a 45.46% reduction in entries. It is not a measured reduction
of that size in memory or total time. No complete global cubic coefficient map
or winning full-family factor was an input to this constructor.

The development sequence was:

| Direct family | Gram entries | Exact local-model interval, mHa |
|---|---:|---:|
| Coherent pair-supported cubics; local widths 2 and 3 | 325,636 | 13.599405179 |
| Also put three-orbital windows in the collective block | 415,556 | 2.482976876 |
| Also add four-orbital local and collective supports | 702,332 | 1.556532063 |

The last interval becomes 1.556539277 mHa after charging rotation error at both
endpoints. The preceding levels are failed target searches, not family-limit
theorems. The final level used an exact-checked operator embedding of the
preceding numerical proposal. Its 118-second preparation, 454-second solve and
240-second local replay do not include all prior work. A further 245-second
complete original-model replay was a charged validation/development expense.
The cold runner avoids this duplicate replay.

The successful development lower used the previously discovered local MPS.
Discovery of that state and every failed preceding proof belong in end-to-end
accounting. The old fully cubic nonsinglet screen was replaced with a directly
constructed quadratic magnetic-sector screen before this final result.

## A cheaper nonsinglet component

The H10 quadratic M_S=1 construction used **88,500 Gram entries** and took
**57.85 seconds** for preparation, numerical search and exact checking. Its
verified lower endpoint lies roughly **15.45 mHa above** the retained ground
upper, so it is sufficient for complete fixed-N certification. It requires
neither an MPS moment input nor a cubic dictionary. This success is specific to
the tested input; the frozen procedure permits local/collective cubic extensions
when the quadratic screen is insufficient.

## Preparation improvement with an exact equivalence check

Relabeling and caching bounded-degree spin patterns reduced the same H10
three-window family's fresh preparation from **320.423 to 71.953 seconds**.
Every coefficient map and numerical preparation input compared bit-identically.
This is a **4.45x preparation-stage speedup**. Preparation peak RSS increased
slightly, from about 672 to 697 MB. It is not a complete solver speedup. The
matched reference is given this same optimization.

## A complete small coupling path

On H4, the exact support partition into two orbital pairs was followed through
lambda=0, 1/4, 1/2 and 1. The exact local-model widths were **0.002783053,
0.007518891, 0.009877641 and 0.012779300 mHa**. The fully coupled endpoint became
**0.012779484 mHa on the original input** after rotation transfer.

All particle-number restrictions were global. Hopping and other charge-transfer
terms were scaled with the rest of the interfragment interaction. The four-point
continuation took **65.60 seconds**, with source integrals and orbital construction
additional. The zero-coupling state was fresh; subsequent states were explicit
continuations. This proves a small interacting example, not favorable scaling.

## Singlet-compatible exact dual repair

On the H4 pair family, the repair enforced **971 equations of rank 140**, including
**27 exact trace-null vectors**, and accepted a **1/1000 trace mixture**. The
constructor and independent checker used occupation counting, not full-sector
enumeration. The complete export/repair/check sequence took **5.84 seconds**.

Its exact dual objective bounds attainable physical-H lower certificates by
approximately **-3.667250370900516 Ha**. Relative to the fixed upper it proves an
interval floor of only **0.250415 mHa**. That does **not** rule out the 1.6 mHa
target. The result is a functional exact repair diagnostic, not an obstruction
to a successful compact method.

The coordinate control also confirmed why canonical support/parity masks cannot
be copied into dense local orbitals. The complete positive/ideal cone transforms
equivalently; a native coefficient-L1 residual penalty is basis dependent. Both
the derivation and that limitation are in [MATHEMATICS.md](MATHEMATICS.md).

The H4 complete-family coordinate comparison used the unrestricted complete
cubic positive/ideal spaces rather than the old index-pruned family. The exact
native-coordinate intervals were **0.000894283 mHa in canonical orbitals** and
**0.004750749 mHa in local orbitals**. Preparation took **3.657 and 2.995 seconds**,
the numerical solves **3.834 and 3.514 seconds**, and local replay **1.182 and
1.557 seconds**, respectively. These are conditional stages, not two fresh
integral-to-certificate pipelines or an equal-optimum theorem for the different
native residual penalties. The local proof does not become uniformly cheaper
merely because the state becomes easier to construct.

Exact receipts are in `results/interacting_scaling_20260915/cases/`, the coupling
receipt in `coupling/h4_coupling/results.json`, the bitwise comparison in
`cases/h10_collective3_patterns/identical_maps.json`, and the repaired dual in
`duals/h4_pair/`.
