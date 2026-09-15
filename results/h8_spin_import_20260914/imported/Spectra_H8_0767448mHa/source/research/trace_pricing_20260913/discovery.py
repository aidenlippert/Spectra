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
    diagonal = np.diag(metric)
    if np.any(diagonal <= 0):
        raise np.linalg.LinAlgError('Nonpositive trace-metric diagonal')
    # Invertible diagonal equilibration removes arbitrary coordinate scales
    # before the definite eigensolve. No ridge or rank truncation is added.
    d = 1/np.sqrt(diagonal)
    gram_scaled = ((gram+gram.T)/2)*d[:, None]*d[None, :]
    metric_scaled = ((metric+metric.T)/2)*d[:, None]*d[None, :]
    values, vectors = eigh(gram_scaled, metric_scaled)
    return values, d[:, None]*vectors


def pricing_basis(coefficients, monomial_trace):
    if coefficients.ndim != 2 or coefficients.shape[0] < coefficients.shape[1] or not np.all(np.isfinite(coefficients)):
        raise ValueError('Invalid pricing coefficient frame')
    Q, R = np.linalg.qr(coefficients, mode='reduced')
    singular = np.linalg.svd(R, compute_uv=False)
    if singular[-1] <= singular[0]*1e-12:
        raise ValueError('Rank-deficient pricing frame; no automatic truncation')
    metric = Q.T@monomial_trace@Q; metric = (metric+metric.T)/2
    np.linalg.cholesky(metric)
    return Q, R, metric


class Model(CoefficientModel):
    def __init__(self, data, tail):
        start = time.monotonic(); super().__init__(data, tail)
        before = time.monotonic()
        self.trace_values = np.array([float(seed(w, self.p['modes'], self.p['particles'])) for w in self.rows])
        self.trace_metrics = []; self.monomial_traces = []; self.pricing_bases = []; statistics = []
        for mapping, C in zip(self.price_maps, self.frame_coefficients):
            monomial_trace = (mapping@self.trace_values).reshape(C.shape[0], C.shape[0])
            metric = C.T@monomial_trace@C
            metric = (metric+metric.T)/2
            if not np.all(np.isfinite(metric)):
                raise ValueError('Nonfinite uniform-trace metric')
            np.linalg.cholesky(metric)
            eigenvalues = np.linalg.eigvalsh(metric)
            Q, R, stable_metric = pricing_basis(C, monomial_trace)
            self.trace_metrics.append(metric)
            self.monomial_traces.append(monomial_trace); self.pricing_bases.append((Q, R, stable_metric))
            statistics.append({'dimension': len(metric), 'minimum_eigenvalue': float(eigenvalues[0]),
                'condition_estimate': float(eigenvalues[-1]/eigenvalues[0]),
                'QR_metric_condition_estimate': float(np.linalg.cond(stable_metric)),
                'QR_reconstruction_relative_error': float(np.linalg.norm(C-Q@R)/np.linalg.norm(C))})
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
        for gid, (frame, mapping, C, monomial_trace, basis_data) in enumerate(zip(self.frames, self.price_maps, self.frame_coefficients, self.monomial_traces, self.pricing_bases)):
            Q, R, metric = basis_data
            monomial_gram = (mapping@y).reshape(C.shape[0], C.shape[0])
            gram = Q.T@monomial_gram@Q
            values, stable_vectors = normalized_modes(gram, metric)
            vectors = np.linalg.solve(R, stable_vectors)
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
                ordinary = float(column@monomial_gram@column); trace_size = float(column@monomial_trace@column)
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
