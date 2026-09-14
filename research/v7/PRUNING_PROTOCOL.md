# Development test of residual-budget pruning

Before method acquisition, test one conventional magnitude pruning rule and one
supplied work-weighted headroom candidate on the unchanged 36 development
problems. No evaluation problem is opened. Allocation is fixed at one quarter
of tolerance, spread uniformly over 24 possible polynomial increments.

If a step omits `d_j` from `G c_(j-1)/j`, its residual contribution is exactly
`-j d_j t^(j-1)` and has integrated l1 bound `||d_j||_1 T^j`. The remaining
highest coefficient supplies the truncation tail. The proposer accounts for
both; the unchanged checker recomputes them from the actual polynomial.

Magnitude pruning orders terms by coefficient magnitude. The candidate orders
by magnitude divided by one plus the number of generated terms in that
column. Computing the column counts and sorting are charged. This is a
supplied heuristic testing potential headroom, not an acquired method.
Compare complete construction, solution, witness and checking elapsed time,
refusals and storage against both magnitude pruning and the strengthened
adaptive Taylor/BFS/Arnoldi baselines. Improvements over a coarse order grid
alone cannot pass. Timings from this initial scan are diagnostic; any positive
candidate needs repeated paired evaluation before authorizing method search.
