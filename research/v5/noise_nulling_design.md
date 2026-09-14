# v5 bounded design: learned noise nulling from raw scalar observations

## Claim boundary

The proposed task is a finite linear-Gaussian experiment in dimension (d=3) or (4). A learner receives only scalar readouts (y=u^T x) for chosen unit controls (u), with every readout and control execution charged. It must discover nuisance-noise directions from calibration data, then choose controls that suppress the discovered mode while acquiring a weaker mode and finally a target signal. The intended claim is that two retained, data-learned directions reduce the cost of a withheld acquisition task after discovery costs are amortized. It is not a claim of general physical law discovery.

## Exact finite source

Use rational parameters and a declared Gaussian source

\[
 x = s\,e_1+n_1v_1+n_2v_2+\xi,
 \qquad y=u^Tx,
\]

where (s,n_1,n_2,\xi) are independent zero-mean Gaussian variables with rational variances, (v_1,v_2) are unknown orthonormal directions, and isotropic readout noise has rational variance. The source is finite-dimensional and exactly specified by its mean/covariance; Gaussian samples are generated only by the evaluator. A calibration context can set one known excitation amplitude or covariance component, but that excitation must be a declared physical source and its scalar readouts must be charged.

Use (d=3) for the first implementation. Draw (v_1) from a finite rational unit-vector set, and draw (v_2) from a finite orthogonal completion. To test genuine transfer, vary a continuous nuisance amplitude or readout variance in held-out evaluation while keeping the direction law fixed. Do not expose (v_i), eigenvalue ordering, covariance, or source labels to the learner.

## Sequential task

Stage 1 presents calibration mixtures in which the dominant unknown noise mode has variance (lambda_1gglambda_2). From raw scalar observations at selected controls, the learner estimates the top noise direction (hat v_1) and a confidence cone. It may then choose (u\perp\hat v_1) to suppress the dominant mode while measuring the target or collecting Stage 2 data.

Stage 2 uses the nulled controls to estimate the weaker direction (hat v_2) or its orthogonal complement. Stage 3 presents a withheld target signal whose raw signal-to-noise ratio is poor under unconstrained controls but improves when controls are orthogonal to both retained directions. The target amplitude/sign is estimated from fresh scalar readouts.

The discoveries must be learned from observations, not supplied as PCA labels. The retained artifact includes the estimated direction, confidence radius, calibration covariance estimate, and the control-selection rule. If the confidence cone is too broad, the learner must abstain or use an exploration control.

## A concrete 3D instance

Use rational eigenvalues

\[
\lambda_1=9,
\qquad \lambda_2=1,
\qquad \sigma_\xi^2=1/4,
\qquad \operatorname{Var}(s)=1/16,
\]

with (v_1,v_2) selected from a hidden finite set of rational rotations. A control exactly orthogonal to (v_1) changes the nuisance variance from (9+1+1/4) to (1+1/4); after also nulling (v_2), it becomes (1/4). The exact orthogonality case is a theorem check, while finite-shot estimates use a Davis–Kahan-style or direct concentration radius for the learned direction. The physical implementation must use estimated controls and therefore include leakage (\lambda_1\sin^2\angle(u,\hat v_1)).

Do not select hidden directions so that a single calibration control identifies them. Require at least two noncommuting or nonparallel probe controls, and withhold one probe orientation from training. This prevents the experiment from reducing to answer lookup.

## Finite certificate

For (m) independent scalar Gaussian readings at a fixed control, the sample mean and second moment have exact Gaussian/chi-square concentration bounds. A simple implementation can use a declared sub-exponential bound for (y^2), with all constants checked numerically and a conservative union bound over controls. The certificate must return a confidence cone for each learned direction and propagate it into an upper bound on residual noise variance:

\[
\operatorname{Var}(u^Tx)
\leq \hat\lambda_1\sin^2\theta_1+\hat\lambda_2\sin^2\theta_2+\hat\sigma_\xi^2+\text{confidence margin}.
\]

The target estimator’s risk certificate uses this bound and fresh target data. Training/calibration observations cannot also serve as target verification. A dense covariance eigensolver may validate the estimate, but the reported bound must be derived from the declared sample concentration and angle estimates.

## Required baselines

Evaluate on the same raw scalar observations and control budget:

- frozen controls learned before Stage 1;
- exact retrieval of prior direction/control pairs;
- stateful greedy PCA/nulling using the same covariance estimator;
- stateful Bayesian posterior over the finite direction family;
- a from-scratch PCA/nulling learner rerun after each stage;
- the proposed retained two-direction learner.

The Bayesian/PCA baselines may tie. A tie is compatible with a finite reusable-representation result and rules out an algorithm-superiority claim.

## Positive and null gates

Predeclare a target variance or mean-squared-error threshold, confidence level, calibration budget, and amortization horizon. A positive finite result requires:

1. Stage 1 direction recovery meets its independent confidence cone.
2. Stage 2 uses the Stage 1 artifact to reduce calibration readouts or target variance relative to reset and frozen baselines.
3. Stage 3 transfers to a withheld direction/amplitude/noise context, with no hidden label or free covariance access.
4. Total scalar-readout cost, including failed probes, direction estimation, control synthesis, and verification, amortizes by the declared horizon.
5. A raw-observation Bayesian/PCA learner receives identical controls and can tie without invalidating the representation claim.

A null result is expected if the eigenvalue gap is too small, direction uncertainty causes leakage, or learning both modes costs more than the target savings.

## Leakage and physical assumptions

Do not reveal the covariance matrix, eigenvectors, principal-component ordering, excitation labels, or a full vector observation. Do not count a control pulse as free if its calibration or execution has a physical cost. Each scalar sensor reading is charged, and simultaneous multi-axis readout requires an explicit apparatus model. The evaluator may know the hidden source; the learner and policy-selection code may not.

The smallest honest conclusion is: “A finite raw-observation learner can estimate and reuse nuisance directions to improve later scalar sensing under a declared Gaussian source and concentration certificate.” It does not establish broad scientific compounding unless the retained direction law transfers to genuinely new coupled contexts and survives the matched Bayesian/PCA comparisons.
