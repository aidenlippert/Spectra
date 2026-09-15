"""Use the same optimized numerical solver for new and reference dictionaries."""
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]


def run(case, tag, seconds, restart=None, mu=2., iterations=None):
    sys.path.insert(0, str(ROOT/'.venv-interacting-libs'))
    from research.gpu_acceleration_20260915 import solve as base
    from research.nvidia_followup_20260915.sparse_quotient import SparseQuotient
    base.Quotient = SparseQuotient
    base.run(case, tag, seconds, mu, restart, 'cpu_evd', iterations, 100, None, iterations is None)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    p.add_argument('tag')
    p.add_argument('--seconds', type=float, required=True)
    p.add_argument('--restart', type=Path)
    p.add_argument('--mu', type=float, default=2.)
    p.add_argument('--iterations', type=int)
    a = p.parse_args()
    run(a.case.resolve(), a.tag, a.seconds, a.restart, a.mu, a.iterations)
