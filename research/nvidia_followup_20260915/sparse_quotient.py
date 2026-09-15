"""Sparse QR implementation of the existing numerical ideal elimination."""
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve_triangular


class SparseQuotient:
    def __init__(self, op):
        import sparseqr
        start = time.perf_counter()
        self.op = op
        F = op.free.tocsc()
        ideals = F[:, 1:]
        maximum_column_norm = np.sqrt(np.asarray(ideals.power(2).sum(axis=0))).max(initial=0.)
        tolerance = 1e-11*max(float(maximum_column_norm), 1.)
        U, R, pivots, rank = sparseqr.qr(ideals, tolerance=tolerance, economy=True)
        self.U = U.tocsc()[:, :rank]
        self.UT = self.U.T.tocsr()
        self.R = R.tocsr()[:rank, :rank]
        self.pivots = np.asarray(pivots[:rank], dtype=int)+1
        f0 = F[:, 0].toarray().ravel()
        v = f0-self.U@(self.UT@f0)
        norm = np.linalg.norm(v)
        if not np.isfinite(norm) or norm < 1e-12*max(1., np.linalg.norm(f0)):
            raise ValueError('Energy coordinate is not independent of the ideals')
        self.y0 = v/(norm*norm)
        self.v = v/norm
        self.C = op.AT(self.y0)
        self.offset = float(op.rhs@self.y0)
        self.rhs = self.project(op.rhs)
        defect = 0.
        for left in range(0, ideals.shape[1], 32):
            columns = ideals[:, left:left+32].toarray()
            defect = max(defect, float(np.max(abs(columns-self.U@(self.UT@columns)))))
        normalization_error = float(np.max(abs(F.T@self.y0-np.r_[1., np.zeros(F.shape[1]-1)])))
        self.record = {'implementation': 'SuiteSparseQR_economy', 'ideal_coordinate_count': ideals.shape[1],
            'numerical_ideal_rank': int(rank), 'dependent_coordinates_removed': ideals.shape[1]-int(rank),
            'QR_absolute_rank_tolerance': tolerance, 'ideal_projector_max_defect': defect,
            'dual_normalization_error': normalization_error, 'Q_nonzeros': self.U.nnz,
            'R_nonzeros': self.R.nnz, 'seconds': time.perf_counter()-start}
        if defect > 1e-8 or normalization_error > 1e-8 or not np.isfinite([defect, normalization_error]).all():
            raise ValueError(('Ill-conditioned numerical ideal elimination', self.record))

    def project(self, y):
        return y-self.U@(self.UT@y)-self.v*np.dot(self.v, y)

    def A(self, Q):
        return self.project(self.op.A(Q))

    def AT(self, y):
        return self.op.AT(self.project(y))

    def recover(self, Q):
        residual = self.op.rhs-self.op.A(Q)
        b = float(self.y0@residual)
        values = np.zeros(self.op.free.shape[1])
        values[0] = b
        remaining = residual-self.op.free[:, 0].toarray().ravel()*b
        values[self.pivots] = spsolve_triangular(self.R, self.UT@remaining, lower=False)
        return values
