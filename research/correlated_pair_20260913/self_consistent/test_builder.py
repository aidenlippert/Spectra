import unittest
raise unittest.SkipTest('Rejected preliminary implementation; replaced by test_fixed_guide and test_algebra. See REJECTED_PROTOTYPES.md.')
import json
from pathlib import Path
from .builder import guide_certificate, paired_cubic_certificate, tau_from_perturbation, tau_cross, self_consistent_update
from fractions import Fraction as F
from experiments.marginal_symbolic import mono

def test_guide_exact_acceptance():
    p=Path('results/molecular_collective_20260913/campaign/h4/fixture.json')
    cert,rec=guide_certificate(json.loads(p.read_text()))
    assert rec['lower'] and len(cert['blocks'])==cert['modes']

def test_bad_moments_refused():
    try: guide_certificate({'modes':4,'particles':2,'hamiltonian':{}}, correlated_moments=[])
    except ValueError: pass
    else: raise AssertionError('bad moments accepted')

def test_paired_cubic_exact_acceptance():
    p=Path('results/molecular_collective_20260913/campaign/h4/fixture.json')
    cert,rec=paired_cubic_certificate(json.loads(p.read_text()))
    assert rec['residual_max_degree'] <= 4
    assert len(cert['blocks']) == 2*cert['modes']

def test_tau_cross_and_iteration():
    v=mono(((1,0),(0,0)))
    w={0:F(2)}; t=tau_from_perturbation(v,w)
    assert tau_cross(t,w)==v
    out,h=self_consistent_update(v,w,iterations=1)
    assert h and out
