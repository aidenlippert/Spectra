"""Operator-size-normalized pricing; the existing solver/checker are inherited."""
import time

import numpy as np
from scipy.linalg import eigh

from research.certificate_scaling.commutator_dual_witness import seed
from research.spin_completion_20260913.discovery import Model as CoefficientModel, DIRECTION_DENOMINATOR


def normalized_modes(gram, metric):
    if gram.shape != metric.shape or gram.ndim != 2 or gram.shape[0] != gram.shape[1]:
        raise ValueError('Pricing Gram and trace metric dimensions differ')
    if not np.all(np.isfinite(gram)) or not np.all(np.isfinite(metric)):
        raise ValueError('Nonfinite pricing Gram or trace metric')
    # Definite generalized eigensolve; no ridge and no rank truncation.
    return eigh((gram+gram.T)/2, (metric+metric.T)/2)


class Model(CoefficientModel):
    def __init__(self, data, tail):
        start = time.monotonic(); super().__init__(data, tail)
        before = time.monotonic()
        self.trace_values = np.array([float(seed(w, self.p['modes'], self.p['particles'])) for w in self.rows])
        self.trace_metrics = []; statistics = []
        for mapping, C in zip(self.price_maps, self.frame_coefficients):
            metric = C.T@(mapping@self.trace_values).reshape(C.shape[0], C.shape[0])@C
            metric = (metric+metric.T)/2
            if not np.all(np.isfinite(metric)):
                raise ValueError('Nonfinite uniform-trace metric')
            np.linalg.cholesky(metric)
            eigenvalues = np.linalg.eigvalsh(metric)
            self.trace_metrics.append(metric)
            statistics.append({'dimension': len(metric), 'minimum_eigenvalue': float(eigenvalues[0]),
                'condition_estimate': float(eigenvalues[-1]/eigenvalues[0])})
        self.construction.update({'pricing_rule': 'uniform_fixed_N_trace',
            'trace_metric_construction_seconds': time.monotonic()-before,
            'trace_metric_entries': sum(T.size for T in self.trace_metrics),
            'trace_metric_statistics': statistics, 'construction_seconds': time.monotonic()-start})

    def price(self, dual, span):
        start = time.monotonic()
        rows = [tuple(tuple(letter) for letter in w) for w in dual['rows']]
        if rows != self.rows or len(dual['values']) != len(rows):
            raise ValueError('Dual proposal uses a different coefficient map')
        y = np.asarray(dual['values'], dtype=float)
        if not np.all(np.isfinite(y)):
            raise ValueError('Nonfinite pricing moments')
        proposals = []; diagnostics = []
        for gid, (frame, mapping, C, T) in enumerate(zip(self.frames, self.price_maps, self.frame_coefficients, self.trace_metrics)):
            gram = C.T@(mapping@y).reshape(C.shape[0], C.shape[0])@C
            gram = (gram+gram.T)/2
            values, vectors = normalized_modes(gram, T)
            old = [C@(np.asarray(e['vector'], dtype=float)/DIRECTION_DENOMINATOR) for e in span if e['group'] == gid]
            accepted = []
            for index, value in enumerate(values):
                if value >= -1e-6 or len(accepted) >= 2:
                    break
                vector = vectors[:, index]; norm = np.linalg.norm(vector)
                if not np.isfinite(norm) or norm <= 0:
                    raise ValueError('Invalid generalized pricing direction')
                z = np.rint(vector/norm*DIRECTION_DENOMINATOR).astype(np.int64)
                v = z.astype(float)/DIRECTION_DENOMINATOR; column = C@v
                if old:
                    basis = np.column_stack(old)
                    residual = column-basis@np.linalg.lstsq(basis, column, rcond=1e-12)[0]
                    if np.linalg.norm(residual) < 1e-7*max(1., np.linalg.norm(column)):
                        continue
                if np.linalg.norm(column) < 1e-10:
                    continue
                ordinary = float(v@gram@v); trace_size = float(v@T@v)
                if trace_size <= 0 or ordinary >= 0:
                    continue
                old.append(column)
                accepted.append({'group': gid, 'vector': list(map(int, z)), 'proposed_moment': float(value),
                    'rounded_ordinary_moment': ordinary, 'rounded_trace_size': trace_size,
                    'rounded_normalized_moment': ordinary/trace_size})
            proposals.append(accepted)
            diagnostics.append({'group': gid, 'dimension': len(values),
                'minimum_proposed_normalized_moment': float(values[0]),
                'negative_eigenvalues': int(np.sum(values < -1e-6)), 'new_independent_candidates': len(accepted)})
        return proposals, {'seconds': time.monotonic()-start, 'frames': diagnostics,
            'scope': 'Generalized trace-normalized floating proposals only. Exact energy replay decides every gain.'}
