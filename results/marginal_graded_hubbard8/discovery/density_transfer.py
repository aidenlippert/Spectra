"""Exact matched upper/lower certificates for three U,t,V target models."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,time,hashlib
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_boundary_transfer import compile_block,contract,enclose,replay_target
from experiments.marginal_boundary_unitary import replay_model
BASE=ROOT/'results/marginal_graded_hubbard8';OUT=BASE/'density_transfer'

def main():
    OUT.mkdir(exist_ok=True)
    source=json.loads((BASE/'tiled_eight_upper/certificate.json').read_text())
    cases=[]
    for U,t,V,a,b,lower in [('4','1','1/2','.374704','.177582','-165801/50000'),
        ('4','1','-1/2','.339851','.179145','-283653/100000'),
        ('3','2/3','1/2','.353530','.174082','-21021/10000')]:
        window={'kind':'local_hubbard_block_v1','sites':6,'U':str(5*F(U)/6),'t':t,'V':V,'lower':lower,
            'onsite_profile':[str(F(x)*F(U)/4) for x in source['lower_window']['onsite_profile']],
            'hopping_profile':[str(F(x)*F(t)) for x in source['lower_window']['hopping_profile']],
            'density_profile':[str(F(x)*F(V)) for x in source['lower_window']['hopping_profile']]}
        cases.append({'target':{'U':U,'t':t,'V':V},'a':a,'b':b,'lower_certificate':window})
    inputs={'upper':source['upper'],'hamiltonian':source['hamiltonian'],'targets':cases,'sites':1000000,'precision_bits':160}
    (OUT/'certificate.json').write_text(json.dumps(inputs,indent=2)+'\n')
    start=time.monotonic();targets=[]
    for case in cases:
        target={k:F(v) for k,v in case['target'].items()}
        result=replay_target(inputs['upper'],inputs['hamiltonian'],case['a'],case['b'],125000,case['lower_certificate'],**target)
        compiled=compile_block(inputs['upper'],inputs['hamiltonian'],**target)
        exact=contract(compiled,case['a'],case['b'],3)
        rounded=enclose(compiled,case['a'],case['b'],3,160)
        if not F(rounded['energy_lower'])<=F(exact['energy'])<=F(rounded['energy_upper']):
            raise ValueError('Target small-chain exact expectation outside enclosure')
        linear=replay_model(inputs['upper'],inputs['hamiltonian'],'57277/250000',1000000,operation='linear_filter',**target)
        linear_interval=enclose(compiled,'57277/250000',0,125000,160)
        if not F(linear_interval['energy_lower'])<=F(linear['upper'])<=F(linear_interval['energy_upper']):
            raise ValueError('Target linear transfer disagrees with previous independent CAR contraction')
        fixed_gate=enclose(compiled,'.356039','.177774',125000,160)
        hi=F(result['upper_per_site'])
        if hi>=F(linear['upper_per_site']) or hi>=F(fixed_gate['upper_per_site']):
            raise ValueError('Target gate did not improve both comparison states')
        if F(result['upper_replay']['width_per_site'])>=F(1,10**30):raise ValueError('Target enclosure too wide')
        result.update(exact_24_site=exact,enclosed_24_site=rounded,
            previous_linear_upper=linear,previous_linear_enclosure=linear_interval,
            unretuned_quadratic_upper=fixed_gate,
            improvement_over_linear_per_site=str(F(linear['upper_per_site'])-hi),
            improvement_from_gate_retuning_per_site=str(F(fixed_gate['upper_per_site'])-hi))
        targets.append(result)
        print(json.dumps({'target':result['target'],'lower':float(F(result['lower_per_site'])),
            'upper':float(hi),'width':float(F(result['width_per_site']))}),flush=True)
    sources=['experiments/marginal_boundary_transfer.py','experiments/marginal_boundary_unitary.py',
        'experiments/marginal_local_hubbard_block.py','experiments/marginal_tiled_upper.py',
        'experiments/marginal_symmetry_moments.py','experiments/marginal_determinant_tree.py',
        'experiments/marginal_transfer_verify.py','experiments/marginal_symbolic.py',
        'experiments/marginal_polynomial_sos.py',
        'results/marginal_graded_hubbard8/discovery/singlet_moment_energy.py',
        'results/marginal_graded_hubbard8/discovery/density_transfer.py',
        'results/marginal_graded_hubbard8/density_transfer/certificate.json']
    receipt={'accepted':True,'seconds':time.monotonic()-start,'targets':targets,
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        'scope':'Three specified nearest-neighbor target models each have a fresh physical upper and matching all-Fock translated-window lower. Contact gates were retuned; source H8 state and lower profile shapes were not optimized for targets. No generic molecular or arbitrary Hamiltonian transfer claim.'}
    (OUT/'independent_replay.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('accepted seconds',receipt['seconds'],flush=True)

if __name__=='__main__':main()
