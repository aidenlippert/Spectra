"""Teacher-free observable reduction with a true controlled-dynamics defect bound.

The basis is algebraic. The only state input is the original charge MPS, used
for four initial scalar moments. No trajectory, reduced embedding, eigenvector,
or occupation-sector list is read. A loose bound is a conditional result.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

from experiments.marginal_hunt_car import add, mono, scale
from experiments.marginal_symbolic import decode, encode, product
from research.correlated_pair_20260913.mps_exact import State
from research.rsi_discovery_20260916.ledger import Ledger, digest

FROZEN_FIXTURE = "22905765c6591d88acb76653755c364a78d7ed8024df4893f0400f1ec91b0edf"
FROZEN_STATE = "c0825551e374612ec50b4b1374b08a8ca60a8fac60ea87542f6c1e191236d07a"


def cadd(a, b):
    return add(a[0], b[0]), add(a[1], b[1])


def cscale(a, s):
    return scale(a[0], s), scale(a[1], s)


def cmul(a, b):
    return (add(product(a[0], b[0]), scale(product(a[1], b[1]), -1)),
            add(product(a[0], b[1]), product(a[1], b[0])))


def drift(h, o):
    real, imag = cadd(cmul(h, o), cscale(cmul(o, h), -1))
    return scale(imag, -1), real  # i[H,O]


def flatten(poly):
    return {(part, w): F(c) for part in (0, 1) for w, c in poly[part].items() if c}


def subtract(a, b, coefficient):
    result = dict(a)
    for key, value in b.items():
        result[key] = result.get(key, 0) - coefficient * value
        if not result[key]:
            del result[key]
    return result


def coordinates(target, basis):
    """Exact real linear projection in Gaussian-rational CAR coefficient space."""
    pivots = {}
    for index, polynomial in enumerate(basis):
        v, coefficients = flatten(polynomial), {index: F(1)}
        for pivot, (row, coordinate) in pivots.items():
            coefficient = v.get(pivot, 0)
            v = subtract(v, row, coefficient)
            coefficients = subtract(coefficients, coordinate, coefficient)
        if not v:
            raise ValueError("Observable dictionary is linearly dependent")
        pivot = min(v, key=lambda k: (len(k[1]), k))
        factor = v[pivot]
        pivots[pivot] = ({k: c / factor for k, c in v.items()},
                         {k: c / factor for k, c in coefficients.items()})
    residual, coordinate = flatten(target), {}
    for pivot, (row, coefficients) in pivots.items():
        value = residual.get(pivot, 0)
        residual = subtract(residual, row, value)
        coordinate = subtract(coordinate, coefficients, -value)
    reconstructed = {}
    for index, coefficient in coordinate.items():
        reconstructed = subtract(reconstructed, flatten(basis[index]), -coefficient)
    if subtract(flatten(target), reconstructed, 1) != residual:
        raise AssertionError("Invalid exact dynamics-defect witness")
    return [coordinate.get(i, F()) for i in range(len(basis))], residual


def exp_upper(x, terms=120):
    """Positive Taylor sum plus geometric upper bound on its omitted tail."""
    x = F(x)
    if x < 0 or x >= terms + 2:
        raise ValueError("Rational exponential envelope exceeded")
    term = total = F(1)
    for k in range(1, terms + 1):
        term *= x / k
        total += term
    tail = term * x / (terms + 1) / (1 - x / (terms + 2))
    return total + tail, tail


def norm_inf(matrix):
    return max((sum(map(abs, row), F()) for row in matrix), default=F())


def upper_round(value, denominator=10**12):
    value = F(value)
    return F(-((-value.numerator * denominator) // value.denominator), denominator)


def apply(matrix, vector):
    return [sum((a * b for a, b in zip(row, vector)), F()) for row in matrix]


def advance(matrix, vector, duration, terms=120):
    term, total = list(vector), list(vector)
    for k in range(1, terms + 1):
        term = [v * duration / k for v in apply(matrix, term)]
        total = [a + b for a, b in zip(total, term)]
    factor, tail = exp_upper(norm_inf(matrix) * duration, terms)
    return total, tail * max(map(abs, vector), default=F()), factor


def construct(data, state):
    if digest(data) != FROZEN_FIXTURE or digest(state) != FROZEN_STATE:
        raise ValueError("This frozen control task requires the original H8 Hamiltonian and MPS")
    h = (decode(data["hamiltonian"], data["modes"], 4), {})
    d = ({((1, p), (0, p)): F(c) for p, c in ((6, 1), (7, 1), (10, -1), (11, -1))}, {})
    w = ({((1, p), (0, q)): F(1) for p, q in ((6, 10), (10, 6), (7, 11), (11, 7))}, {})
    y = ({}, {((1, p), (0, q)): F(c) for p, q, c in ((6, 10, -1), (10, 6, 1), (7, 11, -1), (11, 7, 1))})
    basis = [(mono(()), {}), d, w, y]
    phases = [(F(-88687, 1000000), F(1, 2)), (F(1, 2), F(1, 2)),
              (F(1, 2), F(1, 2)), (F(311663, 1000000), F(1, 2))]
    duration = F(1, 2)
    matrices, residuals, witnesses = [], [], []
    for generator in (h, d, w):
        matrix, errors, records = [], [], []
        for observable in basis:
            coordinate, residual = coordinates(drift(generator, observable), basis)
            matrix.append(coordinate)
            # Every CAR monomial has operator norm <=1; |re|+|im| bounds its coefficient magnitude.
            errors.append(sum(map(abs, residual.values()), F()))
            records.append([[part, [list(letter) for letter in word], str(c)]
                            for (part, word), c in sorted(residual.items())])
        matrices.append(matrix); residuals.append(errors); witnesses.append(records)
    physical = State(data, state)
    initial = []
    for real, imag in basis:
        im = physical.expectation(imag) if imag else F()
        if im:
            raise ValueError("A Hermitian observable had a non-real exact expectation")
        initial.append(physical.expectation(real) if real else F())
    predicted, integration_error, defect_error = initial, F(), F()
    phase_records = []
    for u, v in phases:
        matrix = [[a + u * b + v * c for a, b, c in zip(ra, rb, rc)]
                  for ra, rb, rc in zip(*matrices)]
        error = max(a + abs(u) * b + abs(v) * c for a, b, c in zip(*residuals))
        lipschitz = norm_inf(matrix)
        predicted, tail, growth = advance(matrix, predicted, duration)
        rounded = [F(round(x * 10**12), 10**12) for x in predicted]
        rounding = max(abs(x - y) for x, y in zip(predicted, rounded))
        predicted = rounded
        growth, tail = upper_round(growth), upper_round(tail + rounding)
        integration_error = upper_round(growth * integration_error + tail)
        increment = error * ((growth - 1) / lipschitz if lipschitz else duration)
        defect_error = upper_round(growth * defect_error + increment)
        phase_records.append({"duration": str(duration), "u": str(u), "v": str(v),
            "matrix_norm_inf": str(lipschitz), "defect_norm_upper": str(error),
            "integrated_defect_upper": str(defect_error), "integrator_error_upper": str(integration_error)})
    total_time = duration * len(phases)
    # Original uncertainty budgets: both actuator errors <=.001, D and W norm<=2,
    # initial trace distance<=.0005, total integrated phase-flip rate<=.001.
    robustness = 16 * total_time * F(1, 1000) + 4 * F(1, 2000) + 4 * F(1, 1000)
    radius = defect_error + integration_error + robustness
    lo, hi = max(F(-2), predicted[1] - radius), min(F(2), predicted[1] + radius)
    if lo > hi:
        raise ValueError("Reduced enclosure contradicts the physical observable range")
    return {"kind": "exact_observable_defect_reduction_v2", "fixture_sha256": digest(data),
        "state_sha256": digest(state), "basis": [{"real": encode(a), "imaginary": encode(b)} for a, b in basis],
        "matrices": [[[str(c) for c in row] for row in matrix] for matrix in matrices],
        "defect_coefficients": witnesses, "defect_norms": [[str(e) for e in row] for row in residuals],
        "initial_moments": list(map(str, initial)), "phases": phase_records,
        "predicted_D": str(predicted[1]), "D_interval": [str(lo), str(hi)],
        "target": "-3/5", "target_proved": hi <= F(-3, 5), "robustness_allowance": str(robustness),
        "dimensions": {"reduced": len(basis), "full_sector_constructed": 0, "full_trajectories_read": 0},
        "status": "verified_encoded_claim" if hi <= F(-3, 5) else "conditional_proof_obligation",
        "scope": "Original four phases, original amplitude<=1/2, original rational H8 and MPS; no control rescaling",
        "remaining_obligation": "Tighter state-dependent defects or a richer algebraic observable space are required if the target interval fails."}


def replay(data, state, certificate):
    actual = construct(data, state)
    if actual != certificate:
        raise ValueError("Changed target, dynamics, matrix, defect or initial-state certificate")
    return actual


def main():
    source, output = map(Path, sys.argv[1:])
    ledger = Ledger(output)
    inputs = output / "inputs"
    inputs.mkdir()
    receipt = {}
    for name in ("fixture.json", "state.json"):
        raw = (source / name).read_bytes()
        (inputs / name).write_bytes(raw)
        receipt[name] = {"source": str(source / name), "sha256": hashlib.sha256(raw).hexdigest()}
    ledger.append("input_snapshot", files=receipt, teacher_trajectories_imported=False)
    data, state = (json.loads((inputs / name).read_text()) for name in ("fixture.json", "state.json"))
    with ledger.measure("observable_reduction_and_exact_initial_moments"):
        certificate = construct(data, state)
    with ledger.measure("exact_replay_from_original_inputs"):
        replay(data, state, certificate)
    corrupt = json.loads(json.dumps(certificate)); corrupt["target"] = "0"
    with ledger.measure("changed_target_refusal"):
        try:
            replay(data, state, corrupt)
        except ValueError:
            mutation_rejected = True
        else:
            mutation_rejected = False
    (output / "certificate.json").write_text(json.dumps(certificate, indent=2) + "\n")
    summary = {k: certificate[k] for k in ("status", "target_proved", "D_interval", "target", "dimensions", "scope", "remaining_obligation")}
    summary.update(costs=ledger.costs(), audit=ledger.audit(), changed_target_rejected=mutation_rejected)
    (output / "result.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
