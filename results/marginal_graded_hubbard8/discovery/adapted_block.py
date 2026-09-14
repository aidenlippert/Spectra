"""Exact adapted-block upper and restricted polynomial-span limit."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_adapted_block import replay,krylov_lower
from experiments.marginal_hopping_filter import replay as replay_original
from experiments.marginal_local_hubbard_block import replay as replay_local
BASE=ROOT/'results/marginal_graded_hubbard8'
OUT=BASE/'adapted_block'

def main():
    OUT.mkdir(exist_ok=True)
    source=json.loads((BASE/'tiled_eight_upper/certificate.json').read_text())
    # Frozen rational proposal; numerical discovery is not used by this replay.
    inputs={k:source[k] for k in ('upper','hamiltonian','lower_window')}
    inputs.update(eta='0.230522',coefficients=[
        '0.679842623986','-1','0.568002069394','-0.163882003584',
        '0.026374894910','-0.002437363648','0.000127627458',
        '-0.000003505934','0.000000039157'],effective_span_lower='-4.22965897',
        original_eta='57277/250000')
    (OUT/'certificate.json').write_text(json.dumps(inputs,indent=2)+'\n')
    start=time.monotonic()
    baseline=replay_original(inputs['upper'],inputs['hamiltonian'],inputs['original_eta'],1000000)
    constant=replay(inputs['upper'],inputs['hamiltonian'],inputs['original_eta'],[1],1000000)
    if constant['upper']!=baseline['upper']:raise ValueError('Constant polynomial baseline changed')
    adapted=replay(inputs['upper'],inputs['hamiltonian'],inputs['eta'],inputs['coefficients'],1000000)
    if F(adapted['upper'])>=F(baseline['upper']):raise ValueError('No exact energy improvement')
    family=krylov_lower(inputs['upper'],inputs['hamiltonian'],inputs['eta'],8,inputs['effective_span_lower'])
    span_width=(F(adapted['effective_energy'])-F(family['effective_energy_lower']))/8
    if not 0<=span_width<F(3,10**9):raise ValueError('Fixed-eta polynomial-span bracket too wide')
    family['achieved_thermodynamic_density']=str(F(adapted['thermodynamic_objective_per_block'])/8)
    family['width_per_site']=str(span_width)
    lower=replay_local(inputs['lower_window'])
    if lower['sites']!=6 or F(lower['U'])!=F(10,3) or F(lower['t'])!=1:
        raise ValueError('Uniform lower-window target mismatch')
    lo=F(lower['lower'])/5-F(2,1000000)
    hi=F(adapted['upper_per_site'])
    if lo>hi:raise ValueError('Inconsistent ground-energy interval')
    chains=[]
    for N in (8,16,24,64):
        r=replay(inputs['upper'],inputs['hamiltonian'],inputs['eta'],inputs['coefficients'],N)
        chains.append(r)
    sources=['experiments/marginal_adapted_block.py','experiments/marginal_boundary_unitary.py',
        'experiments/marginal_hopping_filter.py','experiments/marginal_tiled_upper.py',
        'experiments/marginal_symmetry_moments.py','experiments/marginal_determinant_tree.py',
        'experiments/marginal_local_hubbard_block.py','experiments/marginal_symbolic.py',
        'experiments/marginal_transfer_verify.py','experiments/marginal_polynomial_sos.py',
        'results/marginal_graded_hubbard8/discovery/singlet_moment_energy.py',
        'results/marginal_graded_hubbard8/discovery/adapted_block.py',
        'results/marginal_graded_hubbard8/adapted_block/certificate.json']
    result={'accepted':True,'seconds':time.monotonic()-start,'adapted':adapted,
        'original_baseline':baseline,'constant_baseline_match':True,'chains':chains,
        'fixed_eta_polynomial_span':family,'lower_replay':lower,
        'million_site_interval':{'lower_per_site':str(lo),'upper_per_site':str(hi),
            'width_per_site':str(hi-lo),'improvement_per_site':str(F(baseline['upper_per_site'])-hi)},
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        'scope':'Physical adapted H8 states with exact full-chain variational energy and fresh six-site lower certificate. Fixed-eta degree<=8 polynomial-span limit only; arbitrary block-state or eta optimization is not certified. Finite H8 computation does not prove general efficient representability.'}
    (OUT/'independent_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'seconds':result['seconds'],'interval':{k:float(F(v)) for k,v in result['million_site_interval'].items()},
                     'fixed_eta_span_width_per_site':float(span_width)}),flush=True)

if __name__=='__main__':main()
