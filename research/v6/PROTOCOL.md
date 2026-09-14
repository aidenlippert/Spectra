# V6 first gate: privileged method headroom

This protocol is written before fitting or scoring the V6 models. V1–V5 code and results remain regression fixtures. No autonomous method acquisition or weight training is authorized by a positive result from a weak current-state baseline alone.

## Proposed method and anticipated obstruction

The diagnostic privilege is a finished procedure: construct a state from finite input/output history, fit a linear update, and compile it into a companion-state realization for prediction under known future forcing. The privilege does not include a true hidden order, coefficients, state trajectory, poles, or the correct model per system.

A strong ARX model can express the same law. The experiment must expose this comparison rather than cap the conventional model below the necessary order. We expect the controlled delay-state representation may tie the conventional model; that is a failed headroom gate for this candidate, even if both beat a current-state predictor. No claim that all possible methods lack headroom follows.

## Families and split

- Synthetic physical family: dissipatively coupled thermal networks, 3–6 hidden temperatures, two controlled heater perturbations, one observed temperature. Public knowledge gives linear time-invariant dynamics, stability, available inputs, noise scale, and order bound 8. The actual graph, conductances, matrices and state values remain evaluator-only.
- Two predeclared noise regimes: white observation noise and temporally correlated observation noise. Independent seeds 62001–62032 generate 32 systems, alternating regimes. Each supplies four training trajectories, two selection trajectories, two interval-calibration trajectories, and two test trajectories. Training/selection/calibration use bounded random step inputs; tests use pulses and sinusoids, with new initial states and inputs. Trajectories have 192 transitions. The sampling unit for cross-system comparisons is the entire system.
- Fixed history orders: 1,2,3,4,6,8. Ridge values: 1e-8,1e-5,1e-2,1e-1. A forecast is 16 steps from an observed eight-sample history. Forecast windows do not overlap in their target outputs. Known future input is shared by all methods; future outputs are unavailable within a forecast.
- Physical anchor: Dowling Lab's independently collected TCLab step and sine records, frozen at commit d250c5b0625afb35007c075df1e0f7125e74017d. Use the first 540 step-record transitions for fitting, next 180 for selection, final 180 for calibration, then the untouched sine experiment for offline forecasts. Only T1 is supplied as output; both heater columns are supplied as inputs. T2 remains unused. Q2 is identically zero, so this is not rich two-input identification or cross-apparatus generalization. Approximately one-second timing jitter is treated as model mismatch, not silently corrected. No active physical experiment is performed.

## Comparisons

1. Privileged compiled history state, with order and regularization selected from observations.
2. Strong stateful ARX using the complete history-order range and the same observations, fitting and selection rules. Its scalar recurrence and the privileged companion realization are implemented separately, with an algebraic equality check.
3. Full-history nearest-neighbor retrieval using every training transition, with validation-selected neighbor count; forecasts retain their working history.
4. Fixed-language sparse history regression, with sequential threshold/refit over the same lag dictionary. This is a restricted linear-library analogue, not a claim to implement all of SINDy or AI Feynman.
5. Current-state ridge model, labelled weak diagnostic.

All arms receive the same raw data. Shared low-level normalization and least-squares routines are permitted; the privilege cannot quietly enter an autonomous arm because no such arm is built at this gate. Established library learning (DreamCoder/Stitch) is deferred until method acquisition is actually justified. N4SID is a relevant later state-space comparison; no placeholder implementation or claim of running it is permitted. An exact matching strong ARX competitor already suffices to disprove strict headroom for this particular compilation method.

## Metrics, costs and gate

Report per-system multi-step normalized RMSE, max-error interval coverage, and validation refusal. The fixed quality target is NRMSE≤0.15 and observed block coverage≥0.90. Intervals are calibrated from held-out blocks using a finite-sample quantile. Different forcing regimes and serial dependence prevent treating those empirical intervals as unconditional distribution-free coverage certificates.

The headroom target is either at least 10% lower repeatable computational work at matched quality/coverage, or strictly better predictive quality/coverage under the identical observation and fit budget, against the strongest stateful comparator. Differences must be paired by independent system. Report bounded-cost concentration bounds only when their conditions hold. Exactly equivalent predictions are an exact zero effect, not a statistically underpowered positive result. A wall-time difference between equivalent implementations alone does not establish discovery headroom; report implementation timing separately.

Count observed input/output records, simulator trajectories and transitions, physical measurements newly performed (zero), model fits/refits, feature/normal-equation work proxies, forecast arithmetic, history storage, retrieval distance work and elapsed fitting/evaluation time separately. Count unsuccessful candidates and validation failures. Do not add simulator samples to hardware measurements or charge historical V4 calibration to V5/V6.

If the privileged method fails this gate, freeze the null and diagnose the bottleneck as candidate-method/family headroom. Do not start Stage A/B/C method acquisition on that basis. A new family or candidate requires a new prospective gate rather than post-hoc tuning of this test set.
