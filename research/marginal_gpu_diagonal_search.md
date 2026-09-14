# GPU-directed compact overlap search

The running Lambda A100 was handed to this task after the user approved
responsible spending up to$2,600 total. It has now been terminated in Lambda;
the dashboard was checked and showed no running instances. The observed
launch-to-termination lifetime cost is approximately$0.72 at the displayed
$1.99/hour rate, including the earlier benchmark period. This is an estimate,
not a provider invoice. See `results/lambda_runs/diagonal_science/lifecycle.json`.

The frozen benchmark matrices were not used as current scientific input.
The current accepted eight-correction v5 certificate was reconstructed
locally from its physical CAR matrices, projector sources and exact diagonal
correction. The tensors cover all94 local sectors, grouped into12 sizes.
Their fresh physical/affine check differed by less than1.7e-15.

There are120 symmetry-invariant reflection-odd diagonal correction shapes.
The current certificate uses eight. The112 unused shapes each fit the
existing64-entry verifier limit when added individually. A small numerical
linear program proposed a direction for each shape while allowing the old
profiles and penalties to adjust. Seven step lengths for each direction,
plus the baseline, produced785 independent candidates.

The A100 evaluated every candidate on all94 sectors in2.8596 seconds after
0.4387 seconds of setup/transfer and1.8974 seconds of warmup.25 independent
scalar CPU checks covered the baseline, winners, a boundary and random
candidates; maximum disagreement was7.501e-13. The GPU evaluated positive
and negative telescoping coefficients in float64, with the correct density
objective `(ell-alpha*theta_half-beta*theta_joint)/5`.

The resident tensors used35,486,336 bytes. The largest assembled matrix batch
was368,640,000 bytes; these counts exclude vendor eigensolver workspace.
No speedup factor is claimed for this scientific run because only a subset
was timed on CPU. The earlier complete1024-candidate benchmark independently
measured about4x versus the same-instance single-thread CPU reference; that
measurement concerns its frozen inputs, not this new search.

The best screened extra shape was basis index31. Its finite screening step
improved periodic density by about5.69e-7. Local CPU SLSQP refinement of that
selected nine-shape family, including fresh CAR reconstruction, produced
-0.6436570569233961 per periodic site in66 matrix evaluations. The fresh standard-library exact energy replay has now accepted all blocks
in401.725 seconds. The million-site open lower is
`-2011436115385613/3125000000000000 = -0.6436595569233962` per site,
with the historical physical upper `-0.6106763470511881`. This uses nine
shapes and60 nonzero entries within the existing cap. The new lower improves
on the eight-shape lower by approximately2.89336e-5 per site. All source
hashes were checked; subsequently modified proof sources are preserved in
`joint_projector/pre_congruence_sources/`. This certifies the candidate,
not the optimum of the expanded family.

Artifacts:

* `results/lambda_runs/diagonal_science/preparation.json` and `science.npz`
* `results/lambda_runs/diagonal_science/screening.json` and `screening.log`
* `results/lambda_runs/diagonal_science/cpu_refinement.log`
* `results/lambda_runs/diagonal_science/exact_candidate.log`
* `results/marginal_graded_hubbard8/joint_projector/signed_density/symmetric_diagonals_9/extra_31/profile_joint_r1_2_certificate.json`

Source/data and downloaded receipt hashes were compared against the remote
files before termination. Both preparation and GPU screening scripts are
saved beside the results. The remote Python/CUDA environment and benchmark
history are documented in `results/lambda_runs/20260912T092845Z/GPU_HANDOFF.md`.
No further messages may be sent to the archived old task; this task owns
future compute decisions. Exact rational acceptance remains on CPU, and no
universal chemistry or accuracy-versus-cost conclusion follows from this run.
