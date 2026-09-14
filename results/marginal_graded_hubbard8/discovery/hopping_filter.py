"""Replay the correlated extensive upper and its exact filter-family limit."""
from pathlib import Path
from fractions import Fraction as F
from math import lcm
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_hopping_filter import replay
from experiments.marginal_local_hubbard_block import replay as replay_local
from experiments.marginal_polynomial_sos import integer_psd

BASE = ROOT/'results/marginal_graded_hubbard8'
OUT = BASE/'hopping_filter'


def main():
    OUT.mkdir(exist_ok=True)
    inputs = json.loads((BASE/'tiled_eight_upper/certificate.json').read_text())
    inputs.update(eta='57277/250000',merge_shift_lower='-28638459/125000000')
    (OUT/'certificate.json').write_text(json.dumps(inputs,indent=2)+'\n')
    start = time.monotonic()
    lower = replay_local(inputs['lower_window'])
    if lower['sites'] != 6 or F(lower['U']) != F(10,3) or F(lower['t']) != 1:
        raise ValueError('The lower window must reproduce bulk U4 t1')
    chains = []
    for N in (8,16,18,24,64,1000000):
        r = N%8
        receipt = replay(inputs['upper'],inputs['hamiltonian'],inputs['eta'],N,
                         inputs['remainders'][str(r)] if r else None)
        lo,hi = F(lower['lower'])/5-F(2,N),F(receipt['upper_per_site'])
        if lo > hi:
            raise ValueError('Inconsistent correlated energy interval')
        receipt.update(lower_per_site=str(lo),width_per_site=str(hi-lo))
        chains.append(receipt)
        print(json.dumps({'sites':N,'lower_density':float(lo),'upper_density':float(hi),
                          'width_density':float(hi-lo)}),flush=True)
    # On span{Phi,-hPhi}, the energy shift matrix is [[0,-1],[-1,delta]].
    # This proves a lower limit for this upper-state family, not for all states.
    delta = F(chains[-1]['boundary_data']['delta'])
    shift_lower = F(inputs['merge_shift_lower'])
    shifted = [[-shift_lower,F(-1)],[F(-1),delta-shift_lower]]
    scale = lcm(*(x.denominator for row in shifted for x in row))
    psd = integer_psd([[int(x*scale) for x in row] for row in shifted])
    shift_upper = F(chains[-1]['merge_energy_shift'])
    if not 0 <= shift_upper-shift_lower < F(1,10**9):
        raise ValueError('Filter-family bracket did not meet the declared precision')
    family = {'accepted':True,'minimum_merge_shift_lower':str(shift_lower),
              'achieved_merge_shift_upper':str(shift_upper),'width':str(shift_upper-shift_lower),
              'matrix_psd':psd,
              'scope':'Exact restricted-family bracket for linear filters I-eta*h on the fixed H8 trial blocks. Independent real eta values may be chosen at different cuts; each shift obeys the same bracket. This is not a global Hubbard ground-energy lower bound.'}
    sources = ['experiments/marginal_hopping_filter.py','experiments/marginal_tiled_upper.py',
        'experiments/marginal_symmetry_moments.py','experiments/marginal_determinant_tree.py',
        'experiments/marginal_local_hubbard_block.py','experiments/marginal_symbolic.py',
        'experiments/marginal_transfer_verify.py','experiments/marginal_polynomial_sos.py',
        'results/marginal_graded_hubbard8/discovery/singlet_moment_energy.py',
        'results/marginal_graded_hubbard8/discovery/hopping_filter.py',
        'results/marginal_graded_hubbard8/hopping_filter/certificate.json']
    result = {'accepted':True,'seconds':time.monotonic()-start,'lower_replay':lower,
        'chains':chains,'filter_family_bracket':family,
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        'scope':'Fresh bounded H8 physical observables, exact charged CAR identities, analytic cluster-merging induction, physical remainder expectation and all-Fock six-site lower replay. No global state, global sector matrix or exponentially large normalization is expanded.'}
    (OUT/'independent_replay.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
