"""Fresh exact upper-witness tiling paired with the six-site lower proof."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_tiled_upper import replay_tiling
from experiments.marginal_local_hubbard_block import replay as replay_local

OUT = ROOT/'results/marginal_graded_hubbard8/tiled_eight_upper'
BASE = ROOT/'results/marginal_graded_hubbard8'


def main():
    OUT.mkdir(exist_ok=True)
    old = json.loads((BASE/'singlet_moment_sharp/certificate.json').read_text())
    upper = {'embedding':{'basis':old['embedding']['basis']},
             'upper_chebyshev_coefficients':old['upper_chebyshev_coefficients']}
    local = json.loads((BASE/'local_energy_chain/certificate.json').read_text())['local']
    remainders = {str(r):{k:local[f'physical{r}'][k] for k in ('sites','U','t','upper_vector')}
                  for r in (2,4,6)}
    inputs = {'upper':upper,'hamiltonian':json.loads((BASE/'hamiltonian.json').read_text()),
              'remainders':remainders,
              'lower_window':json.loads((BASE/'weighted_window_family/six_site_certificate.json').read_text())}
    (OUT/'certificate.json').write_text(json.dumps(inputs,indent=2)+'\n')
    start = time.monotonic()
    lower = replay_local(inputs['lower_window'])
    if lower['sites'] != 6 or F(lower['U']) != F(10,3) or F(lower['t']) != 1:
        raise ValueError('Six-site overlap window must reproduce bulk U4 t1')
    periodic = F(lower['lower'])/5
    receipts = []
    for N in (10,12,14,64,1000000):
        r = N%8
        receipt = replay_tiling(upper,inputs['hamiltonian'],N,remainders[str(r)] if r else None)
        lower_density = periodic-F(2,N)
        upper_density = F(receipt['upper_per_site'])
        if lower_density > upper_density:
            raise ValueError('Inconsistent exact upper/lower interval')
        receipt.update(lower_per_site=str(lower_density),width_per_site=str(upper_density-lower_density))
        receipts.append(receipt)
        print(json.dumps({'sites':N,'lower_density':float(lower_density),
                          'upper_density':float(upper_density),
                          'width_density':float(upper_density-lower_density)}),flush=True)
    sources = ['experiments/marginal_tiled_upper.py','experiments/marginal_symmetry_moments.py',
        'experiments/marginal_determinant_tree.py','experiments/marginal_local_hubbard_block.py',
        'experiments/marginal_symbolic.py','experiments/marginal_transfer_verify.py',
        'experiments/marginal_polynomial_sos.py',
        'results/marginal_graded_hubbard8/discovery/singlet_moment_energy.py',
        'results/marginal_graded_hubbard8/discovery/tiled_eight_upper.py',
        'results/marginal_graded_hubbard8/tiled_eight_upper/certificate.json']
    result = {'accepted':True,'seconds':time.monotonic()-start,
        'lower_replay':lower,'chains':receipts,
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        'scope':'Fresh exact H8 polynomial upper replay for every chain, physical remainder CAR expectations, and all-Fock six-site local lower replay. Fixed block complexity is independent of global N. The finite open-chain lower uses only the overlap window here; older tiling lower bounds may be stronger for short chains.'}
    (OUT/'independent_replay.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
