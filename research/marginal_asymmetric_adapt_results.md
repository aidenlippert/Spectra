# Tightening certificates after breaking flavor permutation symmetry

The asymmetric ten-mode test now has a certified interval of width **3.109081196939883e-5** at epsilon=1/1000. This is about 193 times tighter than the previous analytic transfer bound. At epsilon=1/100 the width is **8.085205765482178e-4**, about 74 times tighter than analytic transfer. These are model energy units, not Hartree.

The construction reuses the compact matched proof's quartic directions, optimizes new cubic squares for the unequal hoppings, and then reweights the resulting positive directions with a linear program. No full asymmetric quartic SDP solution supplies those directions. A many-body eigenvector is used only to propose the separately checked upper witness.

## New Hamiltonians and retained symmetries

The perturbation is `-sum_i epsilon*(i-2)*(a_i† a_(i+5) + h.c.)` for i=0,...,4. All five hopping coefficients are distinct. The method retains individual pair charges and simultaneous left/right exchange. It does not average the new Hamiltonian over flavor permutations.

A positive square averaged over flavor permutations remains positive even when H lacks that symmetry. Thus the old 21 averaged quartic directions can be retained, while new cubic squares are optimized using only the symmetries the new H actually has.

The resulting optimization has 866 coefficient equations, 316 number-multiplier directions, 432 cubic PSD blocks (maximum dimension 18, 5,551 scalar PSD variables), and 21 nonnegative scalar weights for the old orbit directions. Scalar weights are capped at 100, an explicit extra search restriction. A coefficient residual is penalized and its scalar component is fixed to zero. Export projects negative numerical Gram eigenvalues to zero and accounts for the resulting error through exact residual replay.

## Accepted comparison

| Method | Width at epsilon=1/1000 | Width at epsilon=1/100 |
|---|---:|---:|
| Analytic local-square transfer | 0.0059997063803 | 0.059809533185 |
| Cubic adaptation + 21 orbit directions | 0.00012999584179 | 0.00083530634820 |
| LP reweighting and common-denominator export | **0.000031090811969** | **0.00080852057655** |
| One certified-operator constraint + cubic squares | 0.00041979116427 | 0.0040103433265 |

The polished intervals are approximately:

```text
epsilon=1/1000: 3.280973687071434 <= E0 <= 3.281004777883403
epsilon=1/100:  3.280006084111149 <= E0 <= 3.280814604687698
```

The first LP reduces 1,196 directions to **15 orbit directions plus 238 direct cubic squares**. The second uses **16 orbit directions plus 241 direct squares**. The LP selects and reweights existing directions; it does not discover new quartic vectors. Its recorded solve times were 2.56 and 1.72 seconds, excluding construction and exact replay. Common-denominator rounding is included in the final residual bounds.

The cubic adaptation solves returned `optimal_inaccurate`. Their exact intervals remain valid, but their numerical objectives are not accepted as exact optima. At epsilon=1/1000 the unpolished export had residual norm 1.04e-4; LP reweighting lowered that to about 1.27e-5. Thus a material part of the first improvement comes from finding better weights and reducing certificate error. The larger-perturbation interval remains substantially wider, so no uniform accuracy or convergence rate is established.

An earlier exploratory run, before the scalar-weight cap and scalar-residual gauge were introduced, produced width 8.1339e-5. That certificate is preserved and independently replayable, but the main table uses the current, explicitly constrained formulation and its polished output. The zero-residual 21-direction formulation produced a solver error; that is not an infeasibility proof.

## Reusing a certified operator directly

There is also a stable way to reuse the entire source proof as one positive constraint. If its certified lower bound is L0, then `P=H0-L0 I` is positive on the fixed-N sector. Search for

`Hnew = b I + alpha P + cubic_SOS + (Nhat-N) X + R`, with alpha >= 0.

The coefficient search then needs only body-degree-three rows: 316 projected equations and 76 multiplier directions. P has the same low-body coefficients as H0. Its positivity is supported by the previously checked certificate, not by an unproved new oracle.

Writing the source identity as `H0=b0 I+S0+(Nhat-N)X0+R0`, the exported proof uses scalar `b+alpha*(b0-L0)`, squares `alpha*S0+cubic_SOS`, multiplier `X+alpha*X0`, and residual `R+alpha*R0`. Because `b0-L0=||R0||_1`, exact replay supports the intended bound. A nonzero source residual is allowed and explicitly charged.

This version returned `optimal` numerical proposals and valid rational intervals at both perturbations. Its intervals are wider than those from 21 separate orbit directions. It improves numerical conditioning and reuses an existing bound in low-body coefficient space; it does not resolve the missing asymmetric positivity directions.

## Construction cost and independent reference

Initially constructing each orbit column by rational permutation averaging took about 35 seconds. Reusing the established Reynolds coefficient projection reduced observed preparation to about 3.4 seconds. Tests compare the projected columns to exact orbit averages. These individual timings are not an asymptotic or controlled end-to-end benchmark.

A separately controlled full-quartic reference finished within a 120-second process deadline. It used only pair charges and global exchange: 1,122 PSD blocks, maximum dimension 96, and 59,181 scalar PSD variables. Its numerical objective was **above a valid upper bound**, so it cannot establish the asymmetric quartic optimum. Its exact exported interval was wider than the smaller adapted certificate. See [the reference audit](marginal_asymmetric_reference_results.md).

## Reproduction and remaining target

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_asymmetric_adapt --epsilon 1/1000 --penalty
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_asymmetric_adapt --epsilon 1/100 --penalty
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_asymmetric_polish results/marginal_asymmetric_adapt/1_1000_penalty/certificate.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_asymmetric_polish results/marginal_asymmetric_adapt/1_100_penalty/certificate.json
python3 -S -m experiments.marginal_symmetry_transfer --verify results/marginal_asymmetric_adapt/1_1000_penalty/polished/certificate.json
```

Use `--aggregate` without `--penalty` for the certified-operator version. Accepted receipts include independently recomputed integer-vector upper bounds bound to the same rational Hamiltonian.

The next mathematical/computational target is to generate new asymmetric quartic directions from the remaining violated moment matrices, while keeping optimization small and certificate errors controlled. None of these finite tests establishes the general marginal cone's boundary or scalable quantum chemistry.

## Subsequent native precision correction

The historical LP figures above used SciPy’s bundled HiGHS defaults. The production polisher now uses native HiGHS with a smaller internal matrix-entry threshold. The same raw epsilon=1/1000 proof now yields an independently replayed interval of width **1.929574335056606e-5**, using 13 orbit and 237 direct squares. The controlled backend experiment isolates coefficient dropping as the source of the earlier LP-to-exact loss; see [precision measurements](marginal_lp_precision_probe_results.md). New asymmetric quartic column generation is now implemented separately; its search results must be judged by exact exported bounds, not by negative eigenvalues alone.
