# Conventional baseline strengthening, before method acquisition

The first V7 development run used a coarse Taylor order grid and rebuilt
witness residuals. Its 144 rows are preserved as preliminary comparator data,
not evidence of a learned method. It found 24/36 certified cases for full
Taylor, 27/36 for BFS projection, and 24/36 each for residual expansion and
Arnoldi within the common budgets.

A conventional adaptive Taylor implementation now checks every order, retains
computed coefficients, and uses the exact recurrence to propose a residual
witness directly. If `c_(j+1)=G c_j/(j+1)`, all lower residual coefficients
vanish and only `-G c_m t^m` remains. The independent checker still recomputes
and checks every coefficient. These standard facts are supplied to the strong
baseline. Any saving over the earlier coarse grid belongs to the baseline and
must not be called acquisition or scientific novelty.

Run this baseline on exactly the same development problems and budgets, before
any new method is searched. Preserve complete elapsed construction/solution/
proposal/checking time and failed run cost. No held-out case is opened. The
previous protocol and results retain their hashes; this amendment explains the
prospective strengthening rather than silently changing their interpretation.
