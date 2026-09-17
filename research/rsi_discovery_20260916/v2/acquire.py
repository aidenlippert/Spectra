"""Generation-one model acquisition; no task-specific winner is baked into the prompt."""
import json
from pathlib import Path
import sys

from research.rsi_discovery_20260916.ledger import Ledger
from research.rsi_discovery_20260916.v2.model import propose

PROMPT = """Spectra scientific primitive acquisition: design a faster exact CAR product kernel.
The canonical anticommutation relations are a_i a_j^dagger = delta_ij - a_j^dagger a_i,
and equal creation or equal annihilation operators square to zero. Creation operators precede
annihilation operators; within either group modes are ascending, with the fermionic permutation sign.
A word is a tuple of (creation_flag 0/1, nonnegative_mode) pairs. A polynomial is a dict word: integer.
Input payload has left (tuple), right (tuple), cache (mutable dict). Each word has length <=3;
the inputs need not already be normal ordered. Output a dict of canonical words to nonzero integers.
The empty tuple is the identity. Avoid floating point. All modes and coincidence patterns are arbitrary.
You may call oracle(left,right), which returns the exact dict using the current slower Fraction normalizer.
Workload: many different products of cubic operators over 8..40 modes, from collective SOS dictionaries.
Current oracle caches concrete left/right pairs (100000 entries) but repeats the same algebra on relabeled modes.
Existing spin-twirl pattern caching and paired anticommutator elimination are already in the baseline;
this primitive must handle general CAR products, including mixed coincidence patterns.
Find a reusable structural compilation/memoization method. Cache is per run and is initially empty;
include its construction cost. The primitive must produce exact results even on arbitrary noncanonical inputs.
The eventual task is to make new coefficient-map algorithms discoverable, so retain a clean general interface.
"""


def main():
    out = Path(sys.argv[1]).resolve()
    ledger = Ledger(out)
    result = propose(ledger, PROMPT, "acquire-car", timeout=150)
    (out / "proposal.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"proposal": bool(result), "audit": ledger.audit(), "costs": ledger.costs()}))


if __name__ == "__main__":
    main()
