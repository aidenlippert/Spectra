# Policy-cover kernel review

I found no mathematical counterexample in sections 1–4 of the bridge or in the
current deterministic kernel, conditional on its explicitly supplied premises.

`_cover` uses a prefix-free path set and checks Kraft equality
`sum(2**-len(path)) == 1`. For a binary tree this is a complete leaf cover; the
lexicographic adjacent-prefix check is sufficient because descendants of a
prefix are contiguous. `box_for` applies the same longest-side split to every
path, so sibling boxes partition their parent even when axes change by depth.

`_bounds` computes the exact infinity-distance from a sample point to a cell,
then applies the declared Lipschitz modulus. Thus each cell interval is valid
when every supplied sample interval encloses the same response and the modulus
premise holds. The active-cell lower objective bound and point-feasible upper
bound in `certify` match the bridge proof: a feasible point can only lie in an
active cell, and an accepted point is feasible under its point interval.

The refusal paths are also sound: incomplete or overlapping covers are rejected;
an empty active set proves infeasibility only within the covered decoded domain;
absence of a feasible point returns `unresolved`; and a wide interval does not
become a positive certificate. `construct` retains a complete cover after each
split, charges one oracle call per sampled midpoint, and replays through
`certify`.

The main limitation is epistemic rather than an arithmetic defect: `Sample`
contains no provenance, confidence level, reset/selection metadata, or physical
applicability evidence. An interval can therefore be adversarially supplied and
the kernel will still certify it. This is consistent with the module docstring
and bridge sections 1–4, which make simultaneous enclosure and executable-policy
decoding premises. It must not be reported as a sample-origin or physical
validation check. Likewise, `construct` trusts the oracle to honor its requested
radius; the returned `requested_radius` is accounting only, not enforcement.

The bridge's sample-origin uncertainty, confidence allocation, and model-error
separation therefore remain external obligations. I found no code-level defect
that would justify weakening the current arithmetic certificate, but adding
provenance fields would be necessary if this kernel is later presented as a
complete statistical or physical gate.
