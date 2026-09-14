# Local charge products and a joint dynamic-programming certificate

The eight-site open Hubbard chain at U=4, t=1 now has an exact Q-complement lower bound **−3.81**, discovered from the bare Hamiltonian without an imported metric or a complete charge-pattern list. The strongest fresh search uses 68 local log-factor parameters plus three doublon-profile parameters, 50 restricted convex solves and 164 generated charge-pattern witnesses. Independent standard-library replay verifies the rounded rational factors with at most 390 DP states in a layer. The recorded CPU construction takes 11.23 seconds.

This is a model-specific weighted-row certificate. It does not yet provide a ground-energy interval or characterize the full physical marginal cone.

## Why change the representation?

The preceding polynomial search overcame its numerical obstruction and generated an exact gamma=−14 certificate. A new finite rational dual calculation then exposed a limitation of that metric family: at gamma=−10, every normalized reflection-tied D-times-quadratic charge metric has a negative weighted-row numerator somewhere. Twenty nonnegative physical-row weights give the exact upper bound

`−46072194339942948016402906266 / 68608590558888942633463061365`

on its minimum numerator, approximately −0.671522. The standard-library verifier reconstructs all 1,106 charge patterns, feature matrices, hopping transitions and the mean normalization before checking the rational identity. This is a finite family obstruction, not a statement about arbitrary metrics or quantum states.

A second finite diagnostic makes every positive charge-pattern weight independent. With A equal to negative diagonal interaction energy plus positive hopping multiplicities, the physical numerator is `(-A-gamma I)w`. Positive integer right and left witnesses respectively verify feasibility at gamma=−2.87 and impossibility at gamma=−2.86. This pins the unrestricted charge-only weighted-row limit between those values. It identifies substantial room beyond the quadratic family, but its 1,106-entry weight table is not the desired compact representation.

## The representation

Use

\[
v(q)=D(q)\prod_i f_i(q_i)\prod_{0<j-i\le R}f_{ij}(q_i,q_j),
\qquad q_i\in\{-1,0,1\}.
\]

All factor entries are positive rationals. Neutrality makes `D=sum(q_i^2)/2` the doublon count. The metric vanishes exactly on the all-singly-occupied valence space P and is strictly positive on its complement Q.

Discovery parameterizes the logarithm of the factor product with onsite `q_i,q_i²` and pair `q_i^a q_j^b`, a,b in {1,2}. For any one charge row, each hopping penalty is a positive multiple of an exponential of a linear function of these parameters. Therefore maximizing the worst weighted-row lower bound is a convex optimization problem in this family. The solver proposes parameters; only exact replay accepts their rationalized tables.

For a fixed total D and an allowed edge transition, write delta for its change in doublon count and rho for the ratio of factor products after and before the hop. The row lower bound is

\[
UD-\frac{t}{D}\sum_{\text{edge transitions}} \text{rate}\,(D+\delta)\rho.
\]

Each rho depends only on factors touching the two changed sites. For range R, its dependency interval has length at most 2R+2. The DP maximizes the **sum** of these edge penalties, retaining prefix charge, squared charge and the last 2R+1 charges. It evaluates an edge when the last site in its dependency interval arrives, then merges prefixes with identical remaining state. Fixing D first handles the global D factor without separating incompatible transition maxima.

The charge maximum is attainable by a balanced-spin assignment: alternate spins within every singly occupied run. Even runs are balanced. The number of odd runs is even because the total number of singles is even; choose half their phases each way. This permits all single-single hopping processes simultaneously. Other charge-edge multiplicities do not depend on these phases. The resulting row bound therefore covers every physical spin assignment.

## Exact results

| Construction | Range | Log parameters | Certified Q bound | Peak DP states | Largest local window |
|---|---:|---:|---:|---:|---:|
| Full charge-table convex proposal | 1 | 44 | −4.55 | 129 | 4 |
| Full charge-table convex proposal | 2 | 68 | −4.13 | 390 | 6 |
| Full charge-table convex proposal | 3 | 88 | −4.01 | 560 | 8 |
| Bare-H DP-guided cutting planes | 2 | 68 | **−4.13** | 390 | 6 |

The range-three window spans the whole eight-site system, so that result is not evidence of local-window compression at this size. The range-two certificate occupies about 12 KB, in contrast to the earlier multi-megabyte expanded polynomial decompositions. These are different proof representations, not a sparsification of the same polynomial certificate.

The first three numerical discoveries explicitly used all 1,106 charge patterns. The fourth starts at uniform factors and obtains four initial patterns from the DP. Each restricted convex solve is checked by exact DP; new worst patterns become cuts. It reached the requested threshold after 51 rounds, retaining 155 generated patterns, in approximately 14 seconds in the recorded run. No previous PF weights, polynomial metric, or positive-atom directions enter that construction. These timings are individual CPU observations, not a controlled hardware speedup claim.

The verifier enumerates small local charge windows, not a global list of neutral patterns. Its state and window counts remain exponential in interaction range; bounded range is an assumption, not a universal complexity result. Discovery can still require many cuts on other targets.

## A count profile without a wider memory window

Replace D in the metric by a positive profile g(D), with g(0)=0. At fixed D, every hopping coefficient becomes `g(D+delta)` and the final row divides by `g(D)`. The dependency window and DP state remain unchanged. Discovery uses `g(D)=D exp(b_D)`, with b_1=0; the extra log ratios remain linear inside each exponential, preserving convexity.

The accepted fresh range-two profile search starts with g(D)=D and uniform local factors. It uses four initial DP witnesses, imports no history, and certifies **−3.81** after 50 rounds and 164 generated patterns. Its exact worst weighted row is approximately −3.80608446. The certificate is 11,799 bytes. `product_profile_fresh_R2/provenance.json` binds the Hamiltonian, generator hash, feature choice and absence of a resumed history.

The full-table profile diagnostic reaches approximately −3.79578256 numerically, with a separate exact −3.81 certificate. A prior −3.5 target was not reached; it is retained as a failed target rather than called impossible by numerical status alone. Earlier resumed profile metadata had a local-only parameter count; the fresh run with explicit provenance supersedes it as the headline discovery result.

## Size probes

The g(D)=D, range-two pipeline also produces independently replayed certificates on ten and twelve sites. Initial targets −6 and −8 were not reached within 60 rounds. Retaining those generated cuts and solving at the weaker targets gives:

| Sites | Certified Q bound | Local parameters | Generated cuts | Peak DP states | Factor certificate bytes |
|---|---:|---:|---:|---:|---:|
| 10 | −7 | 88 | 246 | 1,130 | 14,953 |
| 12 | −10 | 108 | 312 | 2,104 | 18,262 |

These are weaker feasibility targets, not equal-accuracy scaling results. The twelve-site balanced spin sector has 853,776 determinants and 73,788 nonvalence charge patterns, counted combinatorially; neither complete list is generated by discovery or replay. The bound's tightness deteriorates with size in these probes. That is a remaining representation/accuracy issue even though verification stays compact at fixed range.

## Validation and next gates

`experiments/marginal_joint_charge_dp.py` binds the certificate to the actual uniform open Hubbard Hamiltonian, validates all rational factor tables and bounds, and refuses unsupported ranges, malformed sectors, exhausted budgets or a missed target. Its four focused tests compare the DP with independent signed fermionic actions on four sites for both D and nonlinear count profiles, verify returned worst-pattern witnesses, recover all four known eight-site D-metric row limits, and check refusal paths. Independent `python -S` replay passes for the accepted product certificates. A read-only independent code review found no defect in the original recurrence or physical multiplicities. The full marginal suite passes **382 tests in 270.359 seconds** after the count-profile extension. Both optional-profile and default-profile discovery paths also pass a four-site executable sanity check.

The next gates are a ground-energy witness and Schur-response integration, transfer beyond the uniform open-chain model, and size/range scaling. Existing energy machinery retains at most 32 determinants through its base oracle, while the complete eight-site balanced valence space has 70. Its response basis is also capped at 32 directions. The standalone DP complement result does not bypass either gate.

An exact rank obstruction now proves that merely connecting the new gap to the existing 32-direction response cannot reach the tested lower endpoints. At gamma=−4.13 and each of tau=−4.20, −4.25 and −4.30, a 33-dimensional principal restriction of the scalar Schur matrix is strictly negative definite. Independent standard-library replay rebuilds all 70 valence CAR action columns, checks H_PP=0, forms their rational leakage Gram and verifies exact LDL pivots. It references 350 determinants and takes 70 source actions; it does not diagonalize the 4,900-state sector. The reported total negative inertias 69, 68 and 67 are numerical diagnostics; the certified minimum response rank is 33.

Any response correction of rank at most 32 has a nonzero null vector within that 33-dimensional negative subspace, so it cannot make the Schur matrix nonnegative. This is an obstruction to this response family at these gamma/tau values, not a physical energy lower bound or a universal representation obstruction. A richer response or a different retained-subspace geometry is necessary. A 14-dimensional singlet valence basis would be a vector-valued embedding, not 14 determinant states; using it also requires a valid ground-spin argument or explicit bounds on the other spin sectors.

The stronger gamma=−3.81 proof does not eliminate this issue. The 70 valence coordinates split into parity classes of sizes 38 and 32 according to the sum of up-spin site indices. Within the 38-coordinate class, the exact leakage Gram is diagonal with entries at least 2. At tau=−4.2, `(-tau)(gamma-tau)=1.638<2`, so the negative scalar-Schur restriction has 38 strictly positive diagonal entries after sign reversal. Independent replay proves a required response rank of at least **38**, with minimum sign-reversed diagonal 181/195. This short parity proof supersedes the more complicated 33-coordinate diagnostic for the new gap.

The earlier least-squares log-factor fits were poor. Their abandoned smooth-minimax implementation had missing D ratios and an incorrect objective sign; its divergence is not mathematical evidence against product metrics. The accepted convex implementation and independent DP checks supersede that diagnostic.

Artifacts: `results/marginal_graded_hubbard8/product_cuts_R2/`, `product_convex_R1/`, `product_convex_R2/`, `product_convex_R3/`, `quadratic_metric_frontier/exact_*`, and `charge_only_limit/`. Reproducible discovery scripts are `discovery/discover_product_cuts.py` and `discovery/optimize_charge_product.py`.


## Subsequent energy integration

The noncoordinate singlet path is now implemented in a bounded research verifier: exact ground interval -4.3 <= E0 < -4.20384, with14 valence-bond vectors and14 degree-five polynomial response vectors. It uses the verified complete singlet embedding, the existing -3.81 DP gap, generalized Gram Schur matrices and an original physical upper witness. The ground-spin step explicitly invokes Lieb's theorem with the original model hypotheses checked. Replay still references4234 of4900 balanced-sector determinants; moment contraction is the next computational bottleneck. See [the energy report](marginal_singlet_energy.md). Earlier response-rank diagnostics remain valid for their specified coordinate reference. The new upper witness rules out -4.2 as a lower energy target.


The subsequent symmetry-orbit moment method tightens the H8 energy interval to [-4.23585,-4.23580593] while retaining the fresh range-two gap. An enumerated range-three count-profile certificate accepts gamma=-3.76 and gives a slightly sharper energy variant. Exact larger-chain moment probes also expose the surviving action-budget limit; see [the current moment report](marginal_symmetry_moments.md).
