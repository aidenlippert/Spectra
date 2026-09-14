# Compact residual certificates and native basis continuation

The residual proof now stores one positive scaling integer per complement coordinate instead of explicitly listing every derived pair factor. The refined H6 certificate is **2,055,818 bytes**, down from **4,592,165 bytes**, and certifies a slightly stronger exact common bound: **−6.530499343207549 Ha**. The required reference threshold is still **−6.26489910104 Ha**.

## One vector encodes the pair factors

The certificate retains the original sparse atoms and diagonal metric W:

\[
R=WH_QW-\sum_a\lambda_a v_av_a^T.
\]

It supplies a positive integer vector t and checks, in exact rational arithmetic,

\[
\gamma\le\min_i\frac{R_{ii}-\sum_{j\ne i}|R_{ij}|t_j/t_i}{W_{ii}^2}.
\]

This is a complete positivity proof. Each residual edge has the implicit positive factor

\[
\frac{|R_{ij}|}{t_it_j}
\left(t_j e_i+\operatorname{sgn}(R_{ij})t_i e_j\right)
\left(t_j e_i+\operatorname{sgn}(R_{ij})t_i e_j\right)^T.
\]

After removing those factors, the displayed inequality certifies the remaining diagonal against γW². The original W stays fixed; t only describes the residual decomposition. Positive common rescaling of t leaves every ratio unchanged.

The H6 proof retains7,422 explicit atoms. Its368 scaling integers describe6,174 and4,503 implicit pair factors in the200- and168-state blocks, respectively. The implicit pairs remain in the same support-two positive cone; this compresses their representation rather than changing the physical Hamiltonian or assuming additional positivity.

The new `spin_scaled_residual_complement_v1` format has its own verifier in `experiments/marginal_scaled_residual.py`. It checks complete physical Q coverage, reconstructs every residual entry, requires bounded positive integer scales, and rejects an overstated threshold. Legacy signed/rational-atom replay does not silently accept the new residual rule. The compact replay uses only the Python standard library.

## Exact preservation required a small refinement

The first compact artifact, `results/marginal_h6/scaled_residual_compact`, used the earlier residual scaling unchanged. Its printed floating bound matched the explicit-pair proof, but its exact bound was lower by approximately8.92594e−17 Ha. That artifact is valid at its own exported bound; it does not exactly preserve the stronger prior endpoint.

The refined constructor increases the weakest scaled row using bounded rational steps. It considers a common multiplication of all integer scales followed by an increment of one coordinate. Every candidate is checked against the full exact minimum, and only strict improvements survive. This is bounded refinement, not a claim of an optimal scaling vector. Scales are capped at10¹²; the original explicit atom amplitude limits are unchanged.

The refined artifact, `results/marginal_h6/scaled_residual_refined`, improves the exact common bound over the explicit-pair certificate by approximately2.42250e−10 Ha while remaining less than half its size. Its block bounds are−6.530499343207549 and−6.479582796747338 Ha. Both compact artifacts explicitly report that the requested reference threshold is unmet.

The numerical scaling proposer was extracted from `marginal_residual_polish.py` for reuse. Compact construction can either propose a scaling numerically or import the previously recorded vector and refine it using only exact arithmetic. A worse scaling preserves the initial certificate's valid bound.

## Persisting a solver basis

The adaptive atom constructor now saves an optional native HiGHS basis alongside its complete direction dictionary. Columns are identified by node or pair; rows by normalization, pair cap, or exact rational atom direction. On reconstruction, semantic keys remap the basis even when pair columns and atom rows appear in a different order.

All prior keys must remain present. New pair columns start nonbasic at their lower bounds, and new inequality rows start with basic slacks. The importer checks key uniqueness, status values, source metadata, and the total basic-variable count against the number of rows. The basis is applied after installing the numerical seed. Only a basis marked valid by HiGHS is saved.

This also handles final priced directions that were saved in the dictionary but not loaded into the previous LP: those become new rows, with any required zero-H pair columns added. A native basis is numerical restart state, not evidence of physical positivity. Exact certificate reconstruction and fallback remain the acceptance gates.

The actual H6 seed saves valid native bases for both blocks, plus256 pending directions per block. The resumed model therefore tests basis extension as well as reordering. The control uses the same source proof and the identical initial dictionary with only the native basis omitted. A canonical hash of that shared state and the paired observations are recorded in `results/marginal_h6/native_basis_benchmark.json`.

| Block | Initial directions | First LP with basis, seconds | First LP without basis, seconds | Numerical first-bound difference |
|---|---:|---:|---:|---:|
| 200 | 8867 | 15.1674 | 64.1589 | 2.87e−10 Ha |
| 168 | 9563 | 19.4472 | 62.4833 | 1.58e−10 Ha |

Both first solves report optimal status. These are single paired observations on the same initial models, about4.23× and3.21× faster for the first LP. Subsequent pricing paths and round budgets differ; concurrent workloads and timing variability preclude a general end-to-end speedup claim.

The six-round resumed search exports exact block bounds−6.538076409565311 and−6.395770555458791 Ha. The200-state block reaches its180-second budget; the168-state block completes six rounds in128.779 seconds. The positive certificate remains separate from the larger dictionary and basis state. Results are in `results/marginal_h6/native_basis_seed`, `native_basis_resume`, and `native_basis_control`.

The resumed certificate occupies2,286,811 bytes; its search state occupies8,887,837 bytes. The55% reduction above concerns the residual positivity certificate, not total discovery storage.

Fresh compact rescaling of the resumed proof proposes weaker exact scales for both blocks, so `results/marginal_h6/native_basis_scaled_residual` preserves the original resumed bounds. It is not stronger than the earlier refined compact certificate. The strongest accepted common bound remains−6.530499343207549 Ha.

## Remaining scope

Compression removes the explicit residual pair list, but the verifier still reconstructs the complete200/168 blocks and references400 spin-sector configurations. The complement threshold is not met. Complete support-four pricing, full-FW4 feasibility, and general scaling remain unresolved. These complement certificates do not add a ground-energy interval to the energy/witness ledger.

The final regression suite passed **273 tests in286.988 seconds**. All six new independent standard-library replays match every saved replay field. New tests cover exact mixed-sign factor reconstruction, the fixed original metric, invalid scales and missing Q coverage, singleton rational arithmetic, strict refinement and fallback, semantic basis reordering, pending fill pairs and atom rows, and malformed basis statuses. The accepted energy/witness ledger remains114.
