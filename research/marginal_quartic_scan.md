# Quartic charged-sector scan

This scan adds degree-four operator polynomials in nonzero conserved hopping-charge sectors to the M10 certificate search. The first candidates are density-dressed hopping operators `n_i h_j` and pure four-annihilator/four-creator sectors. Discovery remains coefficient-space SDP; no many-body occupation matrix is constructed.

The density-dressed Clarabel probe failed numerically before producing a certificate, indicating poor conditioning at this unscaled quartic basis. The pure quartet probe completed the process but did not emit a replay artifact under the initial output path; it is therefore not treated as a result. These are engineering outcomes, not evidence that the sectors are mathematically useless.

The next run should use `marginal_adaptive.assemble_and_solve` with charge blocks and scaling, then export with `operator_degree=4` and replay through the exact symbolic verifier. A certificate is counted only after exact replay.

Update: the adaptive charged2 scan (one creation plus three annihilators and adjoints, together with same-charge pair annihilators) completed for M10,N5,t=1/5. Exact replay gives lower 3.26690228623432, upper 3.28100669557401, width 0.01410440933969. This is a substantial improvement over the cubic width (~0.02385), but the remaining gap is still nonzero.
