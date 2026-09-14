# Operator dictionary pricing: bounded scaling experiment

This experiment is a deliberately modest test of the proposed structural
direction.  It uses an open spinful Hubbard chain at half filling (6 and 8
spin orbitals), and a dictionary of one-body hopping operators.  A dictionary
selects a subspace by applying its operators to the full fixed-particle Fock
space.  The adaptive strategy greedily prices every unused operator by the
lowest Rayleigh eigenvalue after adding it; fixed uses dictionary order and
random uses a seeded matched-size sample.

The receipt is [operator pricing receipt](/Users/aidenlippert/Documents/Spectra/results/certificate_scaling/operator_pricing/receipt.json), and the runnable source is [operator_pricing.py](/Users/aidenlippert/Documents/Spectra/research/certificate_scaling/operator_pricing.py).

The current bounded run used budget 6.  At that budget, adaptive reduced the
numerical projected-energy error from 1.376 to 0.361 Ha for the 6-orbital
case, and from 3.278 to 2.381 Ha for the 8-orbital case, relative to the exact
full-space diagonalization.  Random happened to beat fixed on residual in the
6-orbital case, so this is evidence for a useful greedy signal, not evidence
of a universal win.  The held-out Hamiltonian adds a 0.27 last-bond hopping
perturbation; adaptive residuals were 1.191 and 2.203, compared with fixed
2.052 and 2.412.  The complete traces, supports, and wall times are in the
JSON receipt.

All reported energies are numerical upper estimates.  The full-space
eigenvalue is a numerical reference, not a certified lower bound; residual
norms are computed independently as `||H psi - E psi||`.  This experiment
therefore tests discovery and transfer, not the existence of a rational SOS
certificate.  Its main failure mode is visible: the adaptive pricing cost
grows with the number of unused operators, and the dictionary remains a
one-body surrogate.  A meaningful next test should price Gram directions or
operator blocks while retaining an exact residual verifier.

For a GPU batch, the minimal workload is the same fixed arrays already emitted
by the script: Hamiltonian `H` of shape `(D,D)`, candidate matrix `V` of shape
`(D,M)`, and a held-out `H_hold` of shape `(D,D)`.  Batched candidate scoring
needs `Q` with shape `(D,k)`, computes orthogonalized candidates
`Z=V-Q(Q.T@V)`, and evaluates the smallest eigenvalue of each projected
`(k+1)x(k+1)` matrix.  For this receipt `D={20,70}`, `M={30,50}`, and `k<=6`;
the same kernel can scale to exported batches without changing the numerical
definition.

## Actual CAR SOS follow-up

[`operator_pricing_sos.py`](/Users/aidenlippert/Documents/Spectra/research/certificate_scaling/operator_pricing_sos.py)
prices existing CAR Gram blocks through the restricted SDP assembler in
`experiments.marginal_coefficient`. It compares matched fixed, seeded-random,
and greedy block sets, then reconstructs the coefficient vector independently
from the returned Gram matrices. Run it with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python -m research.certificate_scaling.operator_pricing_sos --modes 4 --budget 2
```

The smoke receipt is
[`sos_receipt.json`](/Users/aidenlippert/Documents/Spectra/results/certificate_scaling/operator_pricing/sos_receipt.json).
The greedy strategy found a feasible two-block certificate with numerical
objective `-0.1403124138` and independently reconstructed maximum coefficient
residual below `1e-10`; the fixed first two blocks were infeasible under the
same SDP, while the seeded random pair was feasible. This is a useful pricing
failure mode: dictionary ordering changes feasibility. It remains a numerical
SDP experiment: solver infeasibility is not an exact cone-separation theorem,
and the residual is not an exact replay certificate.
