"""Recompute bounded connected moments and exact polynomial Rayleigh witnesses.

Discovery uses SciPy for four coefficients; --replay uses the standard library
and recomputes every moment and exact quotient from the declared local model.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_connected_dimer_moments import compile_moments

OUT = ROOT / 'results/marginal_graded_hubbard8/connected_dimer_moments'


def rayleigh(moments, coefficients):
    coefficients = list(map(F, coefficients))
    norm = sum(a*b*moments[i+j] for i,a in enumerate(coefficients)
               for j,b in enumerate(coefficients))
    energy = sum(a*b*moments[i+j+1] for i,a in enumerate(coefficients)
                 for j,b in enumerate(coefficients))
    if norm <= 0:
        raise ValueError('Polynomial boundary must have positive norm')
    return norm, energy/norm


def main(replay=False):
    if replay:
        artifact = json.loads((OUT/'receipt.json').read_text())
        local_moments, _ = compile_moments(2, 4, 1, 8)
        _, dimer_upper = rayleigh(local_moments, artifact['dimer_product_coefficients'])
        checks = []
        for row in artifact['chains']:
            moments, work = compile_moments(row['sites'], 4, 1, 8)
            if list(map(str,moments)) != row['moments']:
                raise ValueError('Recomputed moments disagree')
            norm, upper = rayleigh(moments, row['polynomial_coefficients'])
            if str(norm) != row['norm'] or str(upper) != row['upper']:
                raise ValueError('Recomputed physical Rayleigh witness disagrees')
            if work != row['work']:
                raise ValueError('Recomputed cluster work disagrees')
            if str(dimer_upper*(row['sites']//2)) != row['dimer_product_upper']:
                raise ValueError('Fixed-charge product energy disagrees')
            checks.append({'sites':row['sites'], 'exact_moments':True,
                           'exact_positive_norm':True, 'exact_upper':True,
                           'exact_dimer_product_upper':True})
        result = {'accepted':True, 'checks':checks,
                  'scope':'Recomputed finite cluster contractions and rational Rayleigh witnesses. Arbitrary chain extension depends on the stated analytic connected-interval lemma; not a machine-checked proof or lower-energy certificate.'}
        (OUT/'stdlib_replay.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result))
        return
    import numpy as np
    from scipy.linalg import eigh
    OUT.mkdir(exist_ok=True)
    rows = []
    for sites in [2,4,6,8,10,12,16,24,64,1000,1000000]:
        start = time.monotonic()
        moments, work = compile_moments(sites,4,1,8)
        # Scale H by an exact integer to keep the numerical proposal conditioned.
        scale = max(1,isqrt(moments[2])+1)
        dimension = 2 if sites == 2 else 4
        G = np.array([[float(F(moments[i+j],scale**(i+j)))
                       for j in range(dimension)] for i in range(dimension)])
        K = np.array([[float(F(moments[i+j+1],scale**(i+j+1)))
                       for j in range(dimension)] for i in range(dimension)])
        values,vectors = eigh(K,G)
        vector = vectors[:,0]/max(abs(vectors[:,0]))
        coefficients = [str(F(int(round(x*10**12)),scale**i))
                        for i,x in enumerate(vector)]
        norm, upper = rayleigh(moments,coefficients)
        row = {'sites':sites, 'moments':list(map(str,moments)), 'work':work,
               'polynomial_coefficients':coefficients,'norm':str(norm),
               'upper':str(upper), 'upper_float':float(upper),
               'upper_per_site':float(upper)/sites,
               'numerical_ritz':float(values[0]*scale),
               'seconds':time.monotonic()-start}
        rows.append(row)
        print(json.dumps({k:row[k] for k in ['sites','upper_float','upper_per_site','seconds']}),flush=True)
    # A product of locally filtered fixed-charge dimers retains extensive
    # energy. Every hopping term between dimers has zero expectation.
    dimer_upper = F(rows[0]['upper'])
    for row in rows:
        product_upper = dimer_upper*(row['sites']//2)
        row['dimer_product_upper'] = str(product_upper)
        row['dimer_product_upper_per_site'] = float(product_upper)/row['sites']
    sources = ['experiments/marginal_connected_dimer_moments.py',
               'experiments/marginal_clifford_moments.py',
               'results/marginal_graded_hubbard8/discovery/connected_dimer_probe.py']
    artifact = {'model':'Uniform open half-filled Hubbard U=4,t=1',
                'boundary':'Normalized product of adjacent valence singlets',
                'dimer_product_coefficients':rows[0]['polynomial_coefficients'],
                'dimer_product_qualification':'Each local polynomialH2 boundary has exactly2particles. Product expectations of inter-dimer hopping vanish, so its total Rayleigh energy is B times the exact local quotient.',
                'chains':rows,
                'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                'scope':'Exact fixed-order scalar moments with analytic connected-interval extension. Explicit rational polynomial states give variational upper bounds only. Numerical Ritz values are discovery diagnostics, not certificates of optimality. No claim that a million-site ground state was solved.'}
    (OUT/'receipt.json').write_text(json.dumps(artifact,indent=2)+'\n')


if __name__ == '__main__':
    main('--replay' in sys.argv)
