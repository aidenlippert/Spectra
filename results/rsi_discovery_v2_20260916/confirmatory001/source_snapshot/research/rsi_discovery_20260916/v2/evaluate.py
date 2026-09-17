"""Isolated evaluation of acquired algebra and generated map algorithms."""
import argparse
import json
from pathlib import Path
import random
import statistics
import time
import traceback

from experiments.marginal_symbolic import word_product
from research.rsi_discovery_20260916.v2.algebra import (oracle, task_words, gram_reference,
    validate_map, map_cost, product_cases, validate_polynomial, ladder_identity)
from research.rsi_discovery_20260916.v2.language import compile_program


def kernel(source, seed=0):
    fn = compile_program(source, {"oracle": oracle})
    cache = {}
    start = time.perf_counter()
    count = 0
    for left, right in product_cases(seed, 600):
        p = fn({"left": left, "right": right, "cache": cache})
        validate_polynomial(p)
        if p != oracle(left, right) or not ladder_identity(left, right, p):
            return {"status": "counterexample_to_encoded_claim", "input": repr((left, right)),
                    "candidate": repr(p), "reference": repr(oracle(left, right))}
        count += 1
    checks = time.perf_counter() - start
    timings = []
    for spatial in (4, 6, 8):
        words = task_words({"spatial": spatial, "family": "collective", "weight": -3})
        expected = gram_reference(words)
        for variant in ("reference", "acquired"):
            samples = []
            for repeat in range(3):
                word_product.cache_clear()
                cache = {}
                normal = oracle if variant == "reference" else lambda l, r: fn({"left": l, "right": r, "cache": cache})
                start = time.perf_counter()
                value = gram_reference(words, normal)
                samples.append(time.perf_counter() - start)
                if value != expected:
                    return {"status": "counterexample_to_encoded_claim", "spatial": spatial,
                            "reason": "Coefficient-map mismatch on an acquisition workload"}
            timings.append({"spatial": spatial, "dictionary_words": len(words), "variant": variant,
                            "cold_seconds": samples, "median_seconds": statistics.median(samples), **map_cost(value)})
    return {"status": "verified_encoded_claim", "claim": "exact CAR products on the recorded tests and coefficient-map inputs",
            "universal_program_correctness_proved": False, "independent_ladder_cases": count,
            "checking_seconds": checks, "benchmarks": timings}


def evaluate(source, specs, primitive=None, repeats=3):
    capabilities = {"oracle": oracle}
    if primitive:
        fn = compile_program(primitive, capabilities)
        capabilities["car_kernel"] = fn
    candidate = compile_program(source, capabilities)
    rows = []
    for spec in specs:
        word_product.cache_clear()
        preparation = time.perf_counter()
        words = task_words(spec)
        preparation = time.perf_counter() - preparation
        start = time.perf_counter()
        reference = gram_reference(words)
        check_build = time.perf_counter() - start
        samples, checking = [], []
        for repeat in range(repeats):
            word_product.cache_clear()
            start = time.perf_counter()
            value = candidate({"words": list(words), "cache": {}})
            samples.append(time.perf_counter() - start)
            start = time.perf_counter()
            validate_map(value, words)
            if value != reference:
                differing = next((w for w in set(reference) | set(value) if reference.get(w) != value.get(w)), None)
                return {"status": "counterexample_to_encoded_claim", "task": spec,
                        "row": repr(differing), "expected": repr(reference.get(differing)),
                        "actual": repr(value.get(differing)), "completed": rows}
            checking.append(time.perf_counter() - start)
        rows.append({"task": spec, "dictionary_words": len(words), "preparation_seconds": preparation,
                     "independent_reference_seconds": check_build, "validation_seconds": checking,
                     "cold_seconds": samples, "median_seconds": statistics.median(samples), **map_cost(reference)})
    return {"status": "verified_encoded_claim", "claim": "exact symmetric Gram coefficient map for every recorded input",
            "universal_program_correctness_proved": False, "tasks": rows}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("request", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    request = json.loads(a.request.read_text())
    start = time.perf_counter()
    try:
        if request["kind"] == "kernel":
            result = kernel(request["source"], request.get("seed", 0))
        else:
            result = evaluate(request["source"], request["tasks"], request.get("primitive"), request.get("repeats", 3))
    except BaseException as error:
        result = {"status": "implementation_failure", "error": type(error).__name__ + ": " + str(error),
                  "traceback": traceback.format_exc()}
    result["complete_worker_seconds"] = time.perf_counter() - start
    a.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
