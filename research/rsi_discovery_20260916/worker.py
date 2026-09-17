"""Candidate execution only: no acceptance decisions or history writes."""
from fractions import Fraction
import json
from pathlib import Path
import statistics
import sys
import time

from research.rsi_discovery_20260916.dialect import compile_candidate


def eigenpair(h, den, subset=False):
    import numpy as np
    from scipy.linalg import eigh
    a = np.array(h, dtype=float) / den
    if subset:
        e, v = eigh(a, subset_by_index=[0, 0], check_finite=True, driver="evr")
    else:
        e, v = np.linalg.eigh(a)
    return float(e[0]), [float(x) for x in v[:, 0]]


def cholesky(h, den, shift, digits):
    import numpy as np
    if type(digits) is not int or not 0 <= digits <= 16:
        raise ValueError("Invalid factor precision")
    factor = np.linalg.cholesky(np.array(h, dtype=float) / den - float(shift) * np.eye(len(h)))
    fd = 10**digits
    return [[int(round(float(factor[i, j]) * fd)) for j in range(i + 1)] for i in range(len(h))]


def integer_gram(factor):
    from flint import fmpz_mat
    n = len(factor)
    if not 0 < n <= 500 or any(len(row) != i + 1 or any(type(v) is not int for v in row) for i, row in enumerate(factor)):
        raise ValueError("Bounded triangular integer factor required")
    matrix = fmpz_mat([row + [0] * (n - len(row)) for row in factor])
    gram = matrix * matrix.transpose()
    return [[int(gram[i, j]) for j in range(n)] for i in range(n)]


def exact_margin(h, den, lower, factor, fd, accelerated=False):
    from math import lcm
    from research.rsi_discovery_20260916.evaluator import reference_gram
    gram = integer_gram(factor) if accelerated else reference_gram(factor)
    d = lcm(den, lower.denominator, fd * fd)
    residual = [[h[i][j] * (d // den) - gram[i][j] * (d // (fd * fd))
                 - (int(lower * d) if i == j else 0) for j in range(len(h))] for i in range(len(h))]
    return Fraction(min(r[i] - sum(abs(v) for j, v in enumerate(r) if i != j)
                        for i, r in enumerate(residual)), d)


def main():
    source_path, request_path, output_path = map(Path, sys.argv[1:])
    request = json.loads(request_path.read_text())
    capabilities = {"rational": Fraction, "eigenpair": eigenpair,
                    "cholesky": cholesky, "integer_gram": integer_gram,
                    "reference_margin": exact_margin}
    if request.get("acquired_methods"):
        capabilities["acquired_margin"] = lambda h, den, lower, factor, fd: exact_margin(h, den, lower, factor, fd, True)
    fn = compile_candidate(source_path.read_text(), capabilities)
    samples = []
    outputs = []
    for payload in request["payloads"]:
        times = []
        previous = None
        for repeat in range(request.get("repeats", 1)):
            # Prevent mutation from changing benchmark inputs on subsequent calls.
            fresh = json.loads(json.dumps(payload))
            start = time.perf_counter()
            output = fn(fresh)
            times.append(time.perf_counter() - start)
            if repeat and output != previous:
                raise ValueError("Nondeterministic output across identical repetitions")
            previous = output
        outputs.append(previous)
        samples.append(times)
    Path(output_path).write_text(json.dumps({"outputs": outputs, "samples_seconds": samples,
                                            "median_seconds": sum(statistics.median(t) for t in samples)},
                                           allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
