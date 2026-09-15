"""A finite noisy mode-acquisition protocol with explicit physical primitives.

The learner sees a public instrument contract and binary scalar-sensor records.
The source's hidden frame is confined to the simulator callback. Guarantees are
conditional on the declared linear source, actuator and readout assumptions.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import exp, sqrt
import numpy as np
from experiments.v5_frame import Q, BASE_LAMBDAS, C, S, frame, dot

ACTUATOR_BOUND = F(1, 10**9)
READOUT_BOUND = F(1, 1000)
PREPARATION_TV = F(1, 10**9)
PREPARATION_TIME = 30 * (1 + 2 * BASE_LAMBDAS[0])
TARGET_ERROR = F(1, 20)
STAGE_ERRORS = (F(1, 399), F(1, 399), F(11, 500))
FAMILY = 'dense-orthogonal-Gaussian-prefix-v1'


@dataclass(frozen=True)
class Contract:
    lineage: tuple
    alpha: F
    actuator_bound: F = ACTUATOR_BOUND
    readout_bound: F = READOUT_BOUND
    preparation_tv: F = PREPARATION_TV
    family: str = FAMILY

    def check(self):
        if (type(self.lineage) is not tuple or not 1 <= len(self.lineage) <= 3
                or any(type(s) is not str or not s for s in self.lineage)):
            raise ValueError('one to three named physical components required')
        if type(self.alpha) is not F or not 1 <= self.alpha <= 2:
            raise ValueError('coupling scale outside proved source class')
        for actual, admitted in ((self.actuator_bound, ACTUATOR_BOUND),
                                  (self.readout_bound, READOUT_BOUND),
                                  (self.preparation_tv, PREPARATION_TV)):
            if type(actual) is not F or not 0 <= actual <= admitted:
                raise ValueError('instrument uncertainty exceeds certified budget')
        if self.family != FAMILY:
            raise ValueError('source family changed')


@dataclass(frozen=True)
class LearnedPrefix:
    lineage: tuple = ()
    signs: tuple = ()
    modes: tuple = ()
    residual: tuple = Q[0]
    failure_bound: F = F(0)
    records: tuple = ()
    vector_coordinate_updates: int = 0


def _record_update(previous, contract, outcome, probe):
    stage = len(contract.lineage)
    sign = 1 if outcome == 0 else -1
    axis = Q[stage]
    mode = tuple(C * x + sign * S * y for x, y in zip(previous.residual, axis))
    residual = tuple(-sign * S * x + C * y for x, y in zip(previous.residual, axis))
    record = {'lineage': contract.lineage, 'alpha': str(contract.alpha),
              'probe': tuple(str(x) for x in probe), 'outcome': outcome,
              'threshold': '20', 'active_prefix': stage,
              'actuator_bound': str(contract.actuator_bound),
              'readout_bound': str(contract.readout_bound),
              'preparation_tv': str(contract.preparation_tv)}
    return LearnedPrefix(contract.lineage, previous.signs + (sign,), previous.modes + (mode,),
                         residual, previous.failure_bound + STAGE_ERRORS[stage-1] + contract.preparation_tv,
                         previous.records + (record,), previous.vector_coordinate_updates + 12)


def acquire_next(read, contract, previous=LearnedPrefix(), error_budget=TARGET_ERROR):
    """One physical probe, chosen using the data-derived residual direction."""
    contract.check()
    if not isinstance(previous, LearnedPrefix) or previous.lineage != contract.lineage[:-1]:
        raise ValueError('learned modes do not match physical component lineage')
    if replay_records(previous.records) != previous:
        raise ValueError('retained representation does not match its evidence')
    if type(error_budget) is not F or not 0 < error_budget < F(1, 2):
        raise ValueError('explicit target error in (0,1/2) required')
    stage = len(contract.lineage)
    predicted_failure = previous.failure_bound + STAGE_ERRORS[stage-1] + contract.preparation_tv
    if predicted_failure > error_budget:
        raise ValueError('one-probe confidence is insufficient; no certified answer')
    probe = tuple(-S * x + C * y for x, y in zip(previous.residual, Q[stage]))
    if dot(probe, probe) != 1:
        raise ValueError('corrupted learned representation')
    outcome = read(stage, probe)
    if type(outcome) is not int or outcome not in (0, 1):
        raise ValueError('binary magnitude-threshold observation required')
    return _record_update(previous, contract, outcome, probe)


def replay_records(records):
    """Strong structured retrieval: reconstruct a prefix from its raw records."""
    result = LearnedPrefix()
    for item in records:
        contract = Contract(tuple(item['lineage']), F(item['alpha']),
                            F(item['actuator_bound']), F(item['readout_bound']),
                            F(item['preparation_tv']))
        contract.check()
        if result.lineage != contract.lineage[:-1] or item['active_prefix'] != len(contract.lineage):
            raise ValueError('record lineage changed')
        expected = tuple(-S * x + C * y for x, y in zip(result.residual, Q[len(contract.lineage)]))
        if tuple(map(F, item['probe'])) != expected or item['threshold'] != '20':
            raise ValueError('record is not a valid nulling experiment')
        if type(item['outcome']) is not int or item['outcome'] not in (0, 1):
            raise ValueError('invalid stored outcome')
        result = _record_update(result, contract, item['outcome'], expected)
    return result


def perturb_probe(probe, radius=ACTUATOR_BOUND):
    """Public fixed actuator perturbation, independent of the hidden frame."""
    x = np.asarray([float(v) for v in probe], dtype=float)
    direction = np.asarray((1., -2., 3., -4.))
    direction -= np.dot(direction, x) * x
    norm = np.linalg.norm(direction)
    if norm:
        x = x + float(radius) / 4 * direction / norm
    return x / np.linalg.norm(x)


def make_source(hidden_signs, contract, seed, *, actual_actuator=None,
                actual_readout=None, stationary=False, allow_contract_violation=False):
    """Evaluator-owned OU source. Returned reader exposes no source parameters.

    Fresh sources start within norm 1 of the origin and relax for a fixed,
    explicitly charged time. Each read has an independent bath realization.
    """
    contract.check()
    if type(seed) is not int or len(hidden_signs) != len(contract.lineage):
        raise ValueError('source dimension/seed mismatch')
    frame(hidden_signs)
    actual_actuator = contract.actuator_bound if actual_actuator is None else actual_actuator
    actual_readout = contract.readout_bound if actual_readout is None else actual_readout
    if (actual_actuator > contract.actuator_bound or actual_readout > contract.readout_bound) and not allow_contract_violation:
        raise ValueError('actual instrument exceeds its public contract')
    streams = {j: np.random.default_rng(seed * 7 + j) for j in range(1, len(hidden_signs)+1)}
    ledger = {'scalar_readouts': 0, 'binary_records': 0, 'source_preparations': 0,
              'source_coordinate_preparations': 0, 'relaxation_time_units': 0,
              'reset_operations': 0, 'read_by_prefix': {}, 'stationary_reference': stationary,
              'intentional_contract_violation': allow_contract_violation}
    records = []
    initial = np.asarray((.25, -.25, .25, -.25))
    def read(prefix, probe):
        if type(prefix) is not int or not 1 <= prefix <= len(hidden_signs):
            raise ValueError('only currently available prefix sources may be read')
        if type(probe) is not tuple or len(probe) != 4 or any(type(x) is not F for x in probe) or dot(probe, probe) != 1:
            raise ValueError('exact finite unit-vector command required')
        if any(max(x.numerator.bit_length(), x.denominator.bit_length()) > 64 for x in probe):
            raise ValueError('probe command precision budget exceeded')
        applied = perturb_probe(probe, actual_actuator)
        modes, residual = frame(hidden_signs[:prefix])
        basis = modes + (residual,) + Q[prefix+1:]
        eigenvalues = [1 + float(contract.alpha) * lam for lam in BASE_LAMBDAS[:prefix]] + [1.] * (4-prefix)
        variance = mean = 0.
        for vector, eigenvalue in zip(basis, eigenvalues):
            u = np.asarray([float(x) for x in vector])
            coefficient = float(np.dot(applied, u))
            decay = 0. if stationary else exp(-PREPARATION_TIME / eigenvalue)
            variance += eigenvalue * (1 - decay * decay) * coefficient * coefficient
            mean += decay * float(np.dot(initial, u)) * coefficient
        raw = mean + sqrt(variance) * float(streams[prefix].normal())
        # A fixed bounded readout shift is included; it is not a source label.
        shifted = raw + float(actual_readout)
        outcome = int(abs(shifted) > 20)
        ledger['scalar_readouts'] += 1; ledger['binary_records'] += 1
        ledger['source_preparations'] += 1; ledger['source_coordinate_preparations'] += 4
        ledger['relaxation_time_units'] += PREPARATION_TIME; ledger['reset_operations'] += 1
        ledger['read_by_prefix'][str(prefix)] = ledger['read_by_prefix'].get(str(prefix), 0) + 1
        records.append({'prefix': prefix, 'probe': tuple(str(x) for x in probe),
                        'outcome': outcome, 'alpha': str(contract.alpha)})
        return outcome
    return read, ledger, records


def resource_lower_bounds():
    cold_error = F(123, 496) - PREPARATION_TV
    expected_cold_reads = (1 - 2 * TARGET_ERROR) * (1 + 2 * cold_error)
    return {'one_read_cold_error_lower_bound': cold_error,
            'adaptive_expected_cold_reads_lower_bound': expected_cold_reads,
            'warm_lineage_error_upper_bound': sum(STAGE_ERRORS) + 3 * PREPARATION_TV,
            'one_acquisition_four_targets_warm_reads': 5,
            'four_cold_targets_expected_reads_lower_bound': 4 * expected_cold_reads,
            'net_expected_read_saving_per_four_target_batch_lower_bound': 4 * expected_cold_reads - 5,
            'relaxation_time_per_read': PREPARATION_TIME}
