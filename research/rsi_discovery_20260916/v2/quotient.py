"""Proof-producing dictionary reduction, including every Hermitian Gram cross term."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys

from experiments.marginal_hunt_car import add, adj, mono, scale
from experiments.marginal_symbolic import canonical, encode, number_shift, product
from research.general_mechanism_20260915.sector_quotient import operator_quotient
from research.interacting_scaling_20260915.dictionary import ideal_basis
from research.rsi_discovery_20260916.ledger import Ledger
from research.rsi_discovery_20260916.v2.ideal import IdealSpan, verify_witness


def reduce_dictionary(operators, modes, particles, span):
    physical = operator_quotient(operators, modes, particles)
    basis, coordinates = [], [[] for _ in operators]
    for group in physical["groups"]:
        for pivot, row in zip(group["pivots"], group["coordinates"]):
            basis.append(operators[group["indices"][pivot]])
            by_global = dict(zip(group["indices"], row))
            for j in range(len(operators)):
                coordinates[j].append(by_global.get(j, F()))
    reduced = [add(*(scale(b, c) for b, c in zip(basis, row))) for row in coordinates]
    witnesses, obstructions = [], []
    for i in range(len(operators)):
        for j in range(i, len(operators)):
            before = product(canonical(adj(operators[i])), operators[j])
            after = product(canonical(adj(reduced[i])), reduced[j])
            if i != j:
                before = add(before, product(canonical(adj(operators[j])), operators[i]))
                after = add(after, product(canonical(adj(reduced[j])), reduced[i]))
            difference = add(before, scale(after, -1))
            result = span.witness(difference)
            record = {"i": i, "j": j, "difference": encode(difference),
                      "coordinates": {str(k): str(c) for k, c in result["coordinates"].items()},
                      "remainder": encode(result["unresolved_polynomial"])}
            if result["equivalent_in_declared_ideal"]:
                if not verify_witness(difference, span.columns, result["coordinates"]):
                    raise AssertionError("Independent ideal witness reconstruction failed")
                witnesses.append(record)
            else:
                obstructions.append(record)
    accepted = not obstructions
    return {"original_size": len(operators), "reduced_size": len(basis),
        "basis": [encode(p) for p in basis], "coordinates": [[str(c) for c in row] for row in coordinates],
        "physical_singlet_reduction_proved": True, "declared_raw_ideal_map_equivalence_proved": accepted,
        "includes_all_Gram_cross_terms": True, "energy_scalar_used_as_ideal": False,
        "witnesses": witnesses, "obstructions": obstructions,
        "PSD_transport": "Q -> C^T Q C, followed by the checked free-column corrections",
        "many_body_states_enumerated": 0,
        "scope": "Exact equivalence of this dictionary map modulo the supplied raw free columns; not a general equality of different SOS relaxations"}


def main():
    ledger = Ledger(Path(sys.argv[1]))
    reports = []
    with ledger.measure("exact_operator_null_and_full_cross_term_ideal_proofs"):
        for modes, particles in ((6, 2), (8, 4), (10, 4)):
            columns = [("number_" + str(i), product(number_shift(modes, particles), p))
                       for i, p in enumerate(ideal_basis(modes, [], max_body=2))]
            span = IdealSpan(columns)
            a, b = mono(((0, 0),)), mono(((0, 2),))
            operators = [a, b, add(a, product(a, number_shift(modes, particles)))]
            proof = reduce_dictionary(operators, modes, particles, span)
            proof.update(modes=modes, particles=particles, supplied_ideal_columns=len(columns))
            reports.append(proof)
        # Deliberately insufficient original ideal: a physical null must not be silently admitted.
        obstructed = reduce_dictionary(operators, modes, particles,
                                       IdealSpan([("N-n", number_shift(modes, particles))]))
        obstructed["experiment"] = "same physical reduction with a deliberately narrower declared ideal"
        reports.append(obstructed)
    result = {"reports": reports, "costs": ledger.costs(), "audit": ledger.audit()}
    (ledger.root / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"outcomes": [{k: r[k] for k in ("original_size", "reduced_size", "declared_raw_ideal_map_equivalence_proved", "many_body_states_enumerated")} for r in reports],
                      "costs": ledger.costs(), "audit": ledger.audit()}))


if __name__ == "__main__":
    main()
