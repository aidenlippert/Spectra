# Audit of the v2 noisy-transfer summary

I inspected `experiments/v2_run.py`, `research/v2/protocol.json`, and `results/v2/summary.json`. The reported exact ideal risks are internally consistent:

- frozen/none: 7/16;
- A-only and B-only: 3/8 each;
- AB: 1/4;
- ideal complementarity: 1/16.

The compiled B channel check is exact for visibility (v=1/2): the post-processing probability is (r=v/(2+v)=1/5), and the resulting binary channel has error (1/(2+v)=2/5). The code checks both transition probabilities rather than only one average.

## Calibration-to-risk claims

The acquisition code obtains A first, then uses the learned A mask to choose the B reference axis. The B certificate is therefore conditional on the A acquisition result. Its union bound is computed over the six A records and six B records, with per-record deltas (10^{-3}). The resulting reported failure bound is approximately (7.46\times10^{-5}), below the protocol's broad (0.01) family allowance.

On the event that both learned maps satisfy their model-conditional certificates, the exhaustive held-out-context calculation gives AB risk exactly (1/4) for every listed world. The unconditional upper bound (1/4+delta/4) is valid only because this evaluator's risk definition lies in ([1/4,1/2]): on the calibration-failure event the worst possible increase over (1/4) is (1/4). It would be invalid for a generic loss in ([0,1]), where the safe bound is (1/4+3delta/4) (or (1/4+delta) if no range refinement is used). The report should keep this dependence explicit.

The analogous reported complementarity lower bound (1/16-delta/4) is valid for this shared protocol, but requires a pointwise coupling lemma rather than a generic loss-range union bound. For a fixed target, let (c_A,c_Bin\{0,1\}) indicate whether the inferred A and B axes are correct. Under uniformly random unknown axes and the stated visibility (v), direct averaging gives
\[
R_0=\tfrac12-\tfrac v8,\quad R_A=\tfrac12-\tfrac{v c_A}{4},\quad R_B=\tfrac12-\tfrac{v c_B}{4},\quad R_{AB}=\tfrac12-\tfrac{v c_Ac_B}{2}.
\]
Therefore (R_A,R_B\ge\tfrac12-v/4) and (R_{AB}\le\tfrac12-vc_Ac_B/2), yielding the conservative bound (J\ge v/8-(v/2)\delta) when the acquisition failure probability is at most δ. For (v=1/2), this is (1/16-delta/4). The exact table is sharper: (J=v/8) when (c_A=c_B) (both correct or both wrong), and (J=-v/8) when exactly one map is wrong. Thus (J\ge v/8-(v/4)\delta) if the probability of a disagreement event is bounded by δ. This protocol-specific proof must be cited; the generic four-risk argument alone is insufficient.

The result also explains why positive J is not an acquisition certificate: both maps can be wrong and still give positive J. The acquisition and map-accuracy gates remain necessary.

## Monte Carlo intervals

`evaluate` computes a Hoeffding radius with `comparisons = 7 methods x 3 worlds = 21` and `family_failure = 0.01`. Conditional on the fixed acquired maps, evaluation draws are independent, so the union bound gives simultaneous coverage for the 21 displayed method/world intervals. Because the calibration maps are frozen before evaluation, integrating this conditional guarantee over calibration randomness preserves coverage. The intervals are procedure-coverage statements, not posterior probabilities that the displayed error is the true risk.

The exact held-out enumeration is stronger than the Monte Carlo interval for this finite protocol: it checks all (57^2=3249) held-out context pairs for each acquired map. The Monte Carlo rows should therefore be labeled as validation of the simulator's sampling implementation, while the exact conditional-risk column carries the finite-context result. Shared random numbers pair methods and are useful for contrasts, but do not create independent observations across methods; the 21-way union accounting correctly avoids assuming such independence.

## Misspecification limit

All guarantees are conditional on the declared GF(2) parity family, channel laws, measurement menu, and calibration model. If the physical source has an omitted interaction or context-dependent noise, a high posterior or a successful finite certificate only identifies the best surviving member of the admitted class. A robust statement must enlarge each null to a set of outcome laws (for example, total-variation radius ε); discrimination is impossible when the effective pairwise gap is at most (2ε). If every admitted model is rejected, the correct conclusion is model-class failure, not discovery of the intended mechanism.

Recommended fix: retain the unconditional AB bound with its explicit loss-range premise, replace or qualify the unconditional (J) bound, and add an exact failed-map worst-case calculation or a proof of coherent four-risk coupling before presenting (1/16-delta/4) as certified.
