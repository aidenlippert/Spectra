"""A new research task: direct Gram maps of polynomial operator representations."""
from fractions import Fraction as F
import json
from pathlib import Path
import random
import statistics
import sys
import time
import traceback

from experiments.marginal_hunt_car import add, adj
from experiments.marginal_symbolic import product, word_product
from research.rsi_discovery_20260916.ledger import Ledger
from research.rsi_discovery_20260916.v2.algebra import oracle, task_words, validate_map, map_cost
from research.rsi_discovery_20260916.v2.language import compile_program
from research.rsi_discovery_20260916.v2.library import MethodLibrary


def operators(spec):
    rng = random.Random(spec["seed"])
    pool = task_words(spec["base"])
    pool = rng.sample(pool, min(len(pool), spec["unique_words"]))
    result = []
    for i in range(spec["operators"]):
        selected = rng.sample(pool, min(len(pool), spec["terms"]))
        coefficients = [rng.choice((-7, -3, -1, 1, 2, 5)) for _ in selected]
        if spec.get("large_integers"):
            coefficients = [c * (2**67 + 13 + i) for c in coefficients]
        result.append(dict(zip(selected, coefficients)))
    if spec.get("relations"):
        result[-1] = {w: int(c) for w, c in add(result[0], {w: 2 * c for w, c in result[1].items()}).items()}
    if spec.get("zero_operator"):
        result[-2] = {}
    return result


def reference(polynomials):
    result = {}
    for i, p in enumerate(polynomials):
        for j in range(i, len(polynomials)):
            q = polynomials[j]
            value = product(adj(p), q)
            if i != j:
                value = add(value, product(adj(q), p))
            for word, coefficient in value.items():
                if coefficient:
                    result.setdefault(word, {})[(i, j)] = int(coefficient)
    return result


BASE_SOURCE = '''def propose(payload):
    operators = payload["operators"]
    result = {}
    for i in range(len(operators)):
        for j in range(i, len(operators)):
            value = {}
            for left, a in operators[i].items():
                dl = tuple((1-c,m) for c,m in reversed(left))
                for right, b in operators[j].items():
                    for w, c in oracle(dl, right).items():
                        value[w] = value.get(w, 0) + a*b*int(c)
            if i != j:
                for left, a in operators[j].items():
                    dl = tuple((1-c,m) for c,m in reversed(left))
                    for right, b in operators[i].items():
                        for w, c in oracle(dl, right).items():
                            value[w] = value.get(w, 0) + a*b*int(c)
            for w, c in value.items():
                if c:
                    result.setdefault(w, {})[(i,j)] = c
    return result
'''


def evaluate(request):
    activation = time.perf_counter()
    capabilities = {"oracle": oracle}
    if request.get("retained"):
        info = request["retained"]
        library = MethodLibrary(Ledger(info["ledger"]))
        method = library.activate(info["method_id"], capabilities)
        def monomial_map(payload):
            cache = payload["cache"].setdefault("retained_map_" + info["method_id"], {})
            return method({"words": payload["words"], "cache": cache})
        capabilities["monomial_map"] = monomial_map
    candidate = compile_program(request["source"], capabilities)
    activation = time.perf_counter() - activation
    records = []
    for spec in request["tasks"]:
        preparation = time.perf_counter()
        polys = operators(spec)
        preparation = time.perf_counter() - preparation
        word_product.cache_clear()
        start = time.perf_counter()
        expected = reference(polys)
        checking = time.perf_counter() - start
        samples, validations = [], []
        for repeat in range(request.get("repeats", 3)):
            word_product.cache_clear()
            start = time.perf_counter()
            actual = candidate({"operators": [dict(p) for p in polys], "cache": {}})
            samples.append(time.perf_counter() - start)
            start = time.perf_counter()
            validate_map(actual, [None] * len(polys))
            if actual != expected:
                row = next(w for w in set(actual) | set(expected) if actual.get(w) != expected.get(w))
                return {"status": "counterexample_to_encoded_claim", "task": spec, "row": repr(row),
                    "expected": repr(expected.get(row)), "actual": repr(actual.get(row)), "completed": records}
            validations.append(time.perf_counter() - start)
        records.append({"task": spec, "operator_count": len(polys),
            "distinct_words": len(set(w for p in polys for w in p)), "preparation_seconds": preparation,
            "independent_reference_seconds": checking, "validation_seconds": validations,
            "cold_seconds": samples, "median_seconds": statistics.median(samples), **map_cost(expected)})
    return {"status": "verified_encoded_claim", "claim": "exact Gram map of every supplied polynomial operator, including all cross terms",
            "tasks": records, "activation_seconds": activation, "universal_program_correctness_proved": False}


def prompt(tasks, method_manifest, action, previous=None, feedback=None):
    text = """Discover an exact algorithm for a NEW Spectra representation task: construct the symmetric
Gram coefficient map directly for a dictionary of polynomial operators, each a sparse integer linear
combination of CAR monomials. This is needed when quotient, spin-adapted or elimination coordinates
replace the original monomial dictionary. Preserve all coefficients; no energy/sector simplification is allowed.
Input payload['operators'] is a list of dicts {word: integer coefficient}; each word is a tuple of
(creation_flag, nonnegative_mode). Word degree<=3. All nonzero monomials across the input have the
same particle-number charge. Empty polynomials are allowed. Input order is arbitrary, coefficients
can exceed floating-point integer precision, and linear dependencies are allowed. payload['cache'] is empty.
For coordinate (i,i), output the normal polynomial O_i^dagger O_i. For (i,j), i<j, output
O_i^dagger O_j + O_j^dagger O_i. Output dict canonical_word -> nonempty dict (i,j): nonzero Python int.
Canonical words have creators first, then annihilators, each group ascending with fermionic signs.
Strip exact zero coefficients. The empty tuple is identity. No floating-point approximations.
The baseline capability oracle(left,right) returns dict canonical_word: rational integer for two words.
It caches concrete pairs with a 100000-entry LRU, cleared before each cold benchmark.
You may implement any pure structural algorithm, reuse shared algebra, eliminate duplicate work,
or derive your own compiler. Algorithms must work beyond the supplied development shapes.
The available method library below is loaded automatically; an empty list means the frozen initial
library. The same base model, target information and limits apply to every arm. The retained method
is a callable implementation, not an answer table. You may ignore it when an alternative is cheaper.
Optimize total cold construction time; activation, checking, failures and acquisition are also charged.
Keep the source concise enough to complete within the proposal budget. New shapes are withheld until freezing.
"""
    text += "\nDevelopment specifications:\n" + json.dumps(tasks)
    text += "\nApplicable method library:\n" + json.dumps(method_manifest)
    text += "\nSelected live research action: " + action
    text += "\nCurrent program:\n" + (previous or BASE_SOURCE)
    text += "\nDevelopment observations:\n" + json.dumps(feedback or [])
    return text


if __name__ == "__main__":
    request, output = map(Path, sys.argv[1:])
    start = time.perf_counter()
    try:
        value = evaluate(json.loads(request.read_text()))
    except Exception as error:
        value = {"status": "implementation_failure", "error": type(error).__name__ + ": " + str(error),
                 "traceback": traceback.format_exc()}
    value["complete_worker_seconds"] = time.perf_counter() - start
    output.write_text(json.dumps(value, indent=2) + "\n")
    print(json.dumps(value))
