# Validation

- Include a focused test in the same change that introduces production behavior.
- Scale verification to risk: run focused tests first, integration checks for changed boundaries, and benchmarks or profiles for performance claims.
- Preserve required validation and refusal paths when simplifying an implementation. A shorter implementation must not weaken a technical gate.
- State what was verified and what was not. Report unsupported or indeterminate outcomes as such; they cannot count as passing results.
