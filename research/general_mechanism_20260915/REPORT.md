# General mechanism: exact operator-action quotient

The broad quantum many-body objective is preserved explicitly in the
[mission](../MISSION.md). The next system size is not an automatic research
objective. The current general obstruction remains economical discovery of a
representation and its terminal positivity proof, with controlled growth under
further elimination and other physical queries.

This pass produced a reusable algebraic component and tested it on the actual
retained dictionaries. It did not solve the general many-body problem, obtain a
new molecular interval, or establish a scaling theorem.

## Exact result

For a dictionary O and the faithful singlet trace, the kernel of
G_ij = tau_0(O_i^dagger O_j) is **exactly** the space of combinations that
annihilate every singlet in the supplied fixed-number sector. Exact rational LDL
gives all those combinations and coordinates on the remaining operator space.

Spatial-charge grouping makes many trace products identically zero. The code
constructs the remaining products directly from CAR words, occupation counting,
and the existing exact spin average. It does not enumerate determinants or load
a Hamiltonian, trial MPS, or coefficient-map matrix. It does inherit the supplied
operator dictionary, whose original discovery/preparation cost is not erased.

The theorem applies to any finite dictionary with an available faithful sector
trace. The implementation accepts real rational, spatial-charge-homogeneous
fermionic dictionaries of degree at most three on an even-particle singlet
sector. This stated implementation domain must not be replaced by “arbitrary
many-body operators.”

Applied to the four large frozen H12 blocks, the construction returned:

| Existing block | Dictionary dimension | Exact action rank | Complete nullity | Largest trace block | Trace products evaluated |
|---|---:|---:|---:|---:|---:|
| 21 | 588 | 564 | 24 | 24 | 3,984 |
| 22 | 216 | 204 | 12 | 11 | 876 |
| 23 | 588 | 564 | 24 | 24 | 3,984 |
| 24 | 216 | 204 | 12 | 11 | 876 |

There are **72 independent null directions across these four blocks**. The
earlier hand-constructed result identified 24 directions in blocks 22 and 24;
this pass also resolves the complete physical kernels in blocks 21 and 23.
“Complete” refers to these supplied dictionaries, not the entire operator
algebra or every block in the molecular proof.

The calculation evaluated 9,720 within-charge trace products. An ungrouped
upper-triangle construction for the same dictionaries would contain 393,204
pairs. That is an exact operation-count comparison, not a measured whole-solver
speedup. The singlet sector has dimension 226,512, computed combinatorially;
none of those sector basis states was listed by the constructor.

## Verification and cost

Four focused tests pass. Their independent ladder-action oracle uses complete
two-electron singlet bases on two and three spatial orbitals. It checks every
trace Gram entry, including omitted cross-charge entries, and the action of
every returned null vector. Additional checks cover the particle-number
condition, non-singlet counterexamples, malformed/inexact inputs, and rejection
of non-PSD rational matrices.

A separate arithmetic read of the saved receipt reconstructed every Gram entry
from its positive rational LDL factors and checked every saved kernel vector and
rank/nullity count. This checks saved arithmetic consistency; it is not a second
independent implementation of the physical trace. The independent physical
oracle is the small-system test described above.

The measured constructor function took **15.921 seconds** from the frozen
dictionaries, including its NumPy import and JSON output. Peak process RSS was
**379,928,576 bytes**. The full receipt is **1,616,106 bytes**. Timing excludes
interpreter startup and pre-run imports, original dictionary discovery and
preparation, tests, and development. This is a component measurement, not a
fresh molecular time-to-certificate comparison. An initial test run failed on a
mistyped import name; it was corrected before the passing tests and molecular
dictionary run. No optimization campaign or external compute was run.

The [exact receipt](../../results/general_mechanism_20260915/sector_quotient.json)
contains the Gram matrices, LDL witnesses, kernels, quotient coordinates, and
input hashes. The [run record](../../results/general_mechanism_20260915/sector_quotient.run.json)
states measurement scope and exclusions.

## What this enables, and the remaining gate

This gives a systematic exact action quotient instead of manually guessing
singlet-null channels. It can supply the correct physical support for a
subsequent positivity or dual-repair calculation, and it extends beyond any
particular molecular coefficients or geometry.

It does **not** prove that these identities are already implied by the current
truncated SDP's free ideal columns. Quotienting that SDP as an equivalent
optimization requires an exact ideal-span proof. Adding omitted valid sector
identities instead requires an explicitly strengthened relaxation and an
accepting checker that validates them. Neither step is claimed completed.

No dual obstruction has been repaired, and no accepted energy endpoint changed.
The preserved best H12 interval remains 2.588188 mHa from the earlier campaign.
Exact redundancy removal alone does not establish discovery of the important
nonzero correlations, small induced response operators, or closure under
dynamics. Those remain the general research questions.

The [mathematical derivation](MATHEMATICS.md) states the quotient theorem, the
ideal-equivalence distinction, the conditional dual-repair consequence, and
the unresolved response/terminal-proof construction. The
[primary-source note](SOURCES.md) separates universal correctness, conditional
complexity barriers, and the scope of the existing methods.

## Reproduce

The core and tests use the standard library and existing project CAR routines:

```bash
python -B -m unittest research.general_mechanism_20260915.test_sector_quotient -v
```

The frozen-dictionary loader additionally needs the existing NumPy environment.
Choose a fresh output filename; the runner refuses to overwrite a receipt:

```bash
.venv-correlated/bin/python -B -m research.general_mechanism_20260915.run_sector_quotient \
  results/acceptance_channels_20260915/h12_cached/prepared \
  results/general_mechanism_20260915/sector_quotient_repeat.json
```

This command performs no energy optimization and does not alter existing
certificates or the production accepting checker.
