"""Fresh physical quadratic-filter energy and large-chain interval replay."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,time
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_boundary_transfer import compile_block,contract,enclose
from experiments.marginal_hopping_filter import replay as linear_replay
from experiments.marginal_adapted_block import replay as adapted_replay
from experiments.marginal_local_hubbard_block import replay as lower_replay
BASE=ROOT/'results/marginal_graded_hubbard8';OUT=BASE/'boundary_transfer'

def main():
    OUT.mkdir(exist_ok=True)
    original=json.loads((BASE/'tiled_eight_upper/certificate.json').read_text())
    adapted=json.loads((BASE/'adapted_block/certificate.json').read_text())
    inputs={k:original[k] for k in ('upper','hamiltonian','lower_window')}
    inputs.update(a='0.356039',b='0.177774',precision_bits=160,adapted_comparator=adapted)
    (OUT/'certificate.json').write_text(json.dumps(inputs,indent=2)+'\n')
    start=time.monotonic();compiled=compile_block(inputs['upper'],inputs['hamiltonian'])
    # A different, analytic merging implementation supplies an exact
    # million-site linear-filter reference for the new interval transfer.
    linear=linear_replay(inputs['upper'],inputs['hamiltonian'],'57277/250000',1000000)
    linear_interval=enclose(compiled,'57277/250000',0,125000,160)
    if not F(linear_interval['energy_lower'])<=F(linear['upper'])<=F(linear_interval['energy_upper']):
        raise ValueError('Large-chain linear filter disagrees with independent exact formula')
    exact_cases=[]
    for q in (1,2,3,8):
        exact=contract(compiled,inputs['a'],inputs['b'],q)
        interval=enclose(compiled,inputs['a'],inputs['b'],q,160)
        if not F(interval['energy_lower'])<=F(exact['energy'])<=F(interval['energy_upper']):
            raise ValueError('Small exact energy lies outside rounded enclosure')
        exact_cases.append({'exact':exact,'enclosure':interval})
    million=enclose(compiled,inputs['a'],inputs['b'],125000,160)
    coarse=enclose(compiled,inputs['a'],inputs['b'],125000,96)
    if max(F(million['energy_lower']),F(coarse['energy_lower']))>min(F(million['energy_upper']),F(coarse['energy_upper'])):
        raise ValueError('Precision enclosures do not overlap')
    if F(million['width_per_site'])>=F(1,10**30):raise ValueError('Energy enclosure too wide')
    previous=adapted_replay(adapted['upper'],adapted['hamiltonian'],adapted['eta'],adapted['coefficients'],1000000)
    hi=F(million['upper_per_site'])
    if hi>=F(previous['upper_per_site']):raise ValueError('Quadratic gate failed to improve the adapted linear state')
    lower=lower_replay(inputs['lower_window'])
    if lower['sites']!=6 or F(lower['U'])!=F(10,3) or F(lower['t'])!=1:
        raise ValueError('Uniform lower target mismatch')
    lo=F(lower['lower'])/5-F(2,1000000)
    if lo>hi:raise ValueError('Inconsistent ground-energy interval')
    sources=['experiments/marginal_boundary_transfer.py','experiments/marginal_boundary_unitary.py',
        'experiments/marginal_hopping_filter.py','experiments/marginal_adapted_block.py',
        'experiments/marginal_tiled_upper.py','experiments/marginal_symmetry_moments.py',
        'experiments/marginal_determinant_tree.py','experiments/marginal_local_hubbard_block.py',
        'experiments/marginal_transfer_verify.py','experiments/marginal_symbolic.py',
        'experiments/marginal_polynomial_sos.py',
        'results/marginal_graded_hubbard8/discovery/singlet_moment_energy.py',
        'results/marginal_graded_hubbard8/discovery/boundary_transfer.py',
        'results/marginal_graded_hubbard8/boundary_transfer/certificate.json']
    receipt={'accepted':True,'seconds':time.monotonic()-start,
        'physical_source_replay':compiled['source_upper_replay'],'physical_block_support':compiled['physical_support'],
        'large_linear_exact_reference':linear,'large_linear_enclosure':linear_interval,
        'exact_cases':exact_cases,'million_site_enclosure':million,'precision96_enclosure':coarse,
        'adapted_comparator':previous,'lower_replay':lower,
        'million_site_ground_interval':{'lower_per_site':str(lo),'upper_per_site':str(hi),
            'width_per_site':str(hi-lo),'improvement_over_adapted_per_site':str(F(previous['upper_per_site'])-hi)},
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        'scope':'Fresh exact physical block tensors, dressed cut energy, and rigorously rounded large-chain transfer, paired with an independent all-Fock lower certificate. The two transfer energy endpoints enclose a trial-state expectation; only its upper endpoint is a ground-energy upper bound. Gate optimality and general representability/scaling are not established.'}
    (OUT/'independent_replay.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'seconds':receipt['seconds'],'ground_interval':{k:float(F(v)) for k,v in receipt['million_site_ground_interval'].items()},
        'rounding_width_per_site':float(F(million['width_per_site']))}),flush=True)

if __name__=='__main__':main()
