# Interference branch: grouped weighted occupation blocks

The main experiment groups local transition operators by occupation-conditioned
source and target labels. Each block is bounded with the exact CAR
anticommutator identity, then the nonnegative block matrix is bounded in a
positive diagonal metric. For block bounds `n_ij`, the replayable rule is

    ||B||² <= R(w) C(w),
    R=max_i sum_j n_ij w_j/w_i,
    C=max_j sum_i n_ij w_i/w_j.

The weights are powers of two and all accepted arithmetic is rational. No
determinant graph, Fock basis, or dense molecular matrix is constructed.

On the frozen fixtures, optimizing the seven occupation-label weights gives:

* H6: squared coupling majorant decreases from 2.9836543137395476 to
  2.5479375453610723 Ha² (ratio 0.8539654).
* H8: it decreases from 7.6219966664427105 to 6.479159921872934 Ha²
  (ratio 0.8500607).

The accepted receipts use a bounded 3^7 = 2,187-candidate exponent search
with exponents in {-1,0,1}; they record the fixture hash, selected weights,
and exact rational endpoint. No unbounded weight search is part of replay.

This improves the existing scalar second-sector penalty but remains negative
for the terminal target when combined with the inherited first-sector gap. It
is therefore a quantitative improvement and a precise diagnosis of the next
loss: the metric preserves interference within each local transition block,
while the final row/column aggregation treats distinct occupation channels
independently. It is not a mathematical obstruction to the family.

The coefficientwise CAR fallback remains in `frame_bound.py` for arbitrary
residual words. It gives H6 `-74.341206001704 Ha` and H8
`-138.052076108063 Ha`, so it is only a refusal-capable fallback.

Tests: `test_frame_bound.py` (2 passed), plus the weighted identity check.
Results are in `results/global_response_20260913/interference/`.
