# Partial operator-square development refinement

The full-square candidate failed complete-cost headroom. Before acquisition,
freeze a narrower supplied candidate: attempt squaring only when the current
ordinary bound is at most 1.1 times tolerance and support is at most 40 terms.
Combine the four anticommuting groups with greatest norm bounds; retain the
ordinary sum bound for the remaining groups. The checker verifies exact
partition coverage, each squared-block product and each remaining group.
This reduces the qualifying three-qubit case from 435 to 190 square pair tests.

This is a new development hypothesis chosen after the previous diagnostic;
it is not an untouched evaluation or acquired method. Use the same seven
randomized paired repeats on all 36 development cases and same-checker ablation.
Every additional bound attempt is charged. The first quadratic checker source
is preserved under `research/v7/frozen/v7_quadratic_certificate_v1.py`; the
extension adds a checked sum of disjoint blocks, without modifying the original
power/Bernstein checker. Neither benchmark uses reserved evaluation instances.
