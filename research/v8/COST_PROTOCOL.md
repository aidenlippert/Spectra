# V8 development cost replay (frozen before execution)

Use all 36 V7 development cases, same Hamiltonian, central-Z observable,
gamma in {0,1/5,2}, T in {1/5,1/2}, tolerance 1/1000, degree cap 24,
live support and cache caps 512. No heldout width-5/7 task is generated.

Five repetitions per physical case. Randomize four arms per case/repetition
with seed 8031: V7 adaptive Taylor, conventional fraction-free Taylor,
partition-only V8 cover-schema Taylor, overlapping-cover V8 Taylor. Charge
fresh Generator construction, recurrence and exact exports, every failed norm
probe, witness construction, independent checker, and refusal work. Record
phase times, degree, bounds, source hashes, all rows and failures. Python
process/import startup is shared and outside per-call elapsed; no compilation
or search result is precomputed for an arm. Integer denominator clearing is a
supplied conventional baseline improvement, not an acquired m1.

The two cover arms use byte-identical residual and norm checkers. A cover is
attempted only when all ordinary bounds fail, the best bound is <=1.1 times
tolerance, and support <=40; use 16 coordinate sweeps. This trigger is a
bounded development diagnostic motivated by the existing residual probe, not
an untouched test or learned selection policy. All incurred probe costs count.

Primary descriptive statistic: paired complete-cost ratio at equal validity,
with all 36 status outcomes reported. Whole-suite total includes refusals.
Five repeated timings estimate implementation noise, not 180 independent
physical problems. No confidence or transfer claim follows from this replay.
Any advantage from known denominator clearing strengthens the baseline for a
later acquisition experiment. Any overlapping-cover benefit must survive the
complete calculation, not just a smaller norm or Taylor degree.

Pre-replay correction: the first execution exposed a checker interface cap
that rejected a valid >128-singleton conventional partition in both cover
arms. That execution is preserved under research/v8/invalidated. The cover
checker now permits 512 groups (ordinary partition compatibility), while the
optimizer's family remains bounded at 128 groups and 2048 incidences. No
physics or tolerance changes. Replay all arms under the corrected shared
checker. All other settings and the seed remain fixed.
