"""Use the same optimized numerical solver for new and reference dictionaries."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]


def magnetic_certificate(certificate):
    """Translate only the alpha-sector equality; retain every exact refusal."""
    from fractions import Fraction as F
    from experiments.marginal_symbolic import add, decode, encode
    if F(certificate['casimir_multiplier']) or certificate['spin_ladder_multiplier']:
        raise ValueError('A magnetic-sector proof cannot use singlet-only ideals')
    m = certificate['modes']
    h = decode(certificate['core']['hamiltonian'], m, 4)
    Y = decode(certificate['alpha_multiplier'], m, 2)
    certificate['core']['hamiltonian'] = encode(add(h, Y))
    certificate.update(magnetization=1, singlet=False, spin_twirl=False)
    return certificate


def run(case, tag, seconds, restart=None, mu=2., iterations=None):
    sys.path.insert(0, str(ROOT/'.venv-interacting-libs'))
    from research.gpu_acceleration_20260915 import solve as base
    from research.nvidia_followup_20260915.sparse_quotient import SparseQuotient
    base.Quotient = SparseQuotient
    design = json.loads((case/'design.json').read_text())
    if design.get('magnetization') == 1:
        original_export = base.export
        def magnetic_export(op, grams, x, folder):
            nn, ns = len(op.meta['number_basis']), len(op.meta['spin_basis'])
            if x[1+nn+ns] != 0 or op.meta['ladder_basis']:
                raise ValueError('A magnetic-sector proof cannot use singlet-only ideals')
            result = original_export(op, grams, x, folder)
            path = folder/'certificate.json'
            certificate = magnetic_certificate(json.loads(path.read_text()))
            path.write_text(json.dumps(certificate, separators=(',', ':'))+'\n')
            result['certificate_bytes'] = path.stat().st_size
            return result
        base.export = magnetic_export
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
