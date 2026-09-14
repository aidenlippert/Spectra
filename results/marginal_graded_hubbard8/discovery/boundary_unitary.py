"""Exact contact-rotation limit and transferred target-energy replays."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_boundary_unitary import replay,replay_model,family_bound,model_shift

BASE=ROOT/'results/marginal_graded_hubbard8'
OUT=BASE/'boundary_unitary'


def main():
    OUT.mkdir(exist_ok=True)
    source=json.loads((BASE/'tiled_eight_upper/certificate.json').read_text())
    inputs={k:source[k] for k in ('upper','hamiltonian')}
    inputs.update(x='0.150424',rotation_family_lower='-0.110426865',gram_free='-0.13881',
                  filter_parameter='57277/250000')
    (OUT/'certificate.json').write_text(json.dumps(inputs,indent=2)+'\n')
    start=time.monotonic()
    rotation=replay(inputs['upper'],inputs['hamiltonian'],inputs['x'],1000000)
    data=rotation['edge_reductions']
    family=family_bound(data,inputs['rotation_family_lower'],inputs['gram_free'])
    achieved=F(rotation['gate_energy_shift']);lower=F(family['minimum_shift_lower'])
    if not 0<=achieved-lower<F(1,10**9):raise ValueError('Unitary family bracket too wide')
    family['achieved_shift']=str(achieved);family['width']=str(achieved-lower)
    filtered=model_shift(data,inputs['filter_parameter'],'linear_filter')
    if lower<=filtered:raise ValueError('Rotation family was not excluded as an improvement')
    family['fresh_filter_shift']=str(filtered)
    family['strict_disadvantage_per_cut']=str(lower-filtered)
    targets=[]
    for U,t,V in (('4','1','1/2'),('4','1','-1/2'),('3','2/3','1/2')):
        for operation,parameter in (('unitary',inputs['x']),('linear_filter',inputs['filter_parameter'])):
            for N in (24,1000000):
                result=replay_model(inputs['upper'],inputs['hamiltonian'],parameter,N,
                                    operation=operation,U=U,t=t,V=V)
                targets.append(result)
                print(json.dumps({'target':result['target'],'operation':operation,'sites':N,
                                  'upper_density':float(F(result['upper_per_site']))}),flush=True)
    sources=['experiments/marginal_boundary_unitary.py','experiments/marginal_tiled_upper.py',
        'experiments/marginal_symmetry_moments.py','experiments/marginal_determinant_tree.py',
        'experiments/marginal_local_hubbard_block.py','experiments/marginal_symbolic.py',
        'experiments/marginal_transfer_verify.py','experiments/marginal_polynomial_sos.py',
        'results/marginal_graded_hubbard8/discovery/singlet_moment_energy.py',
        'results/marginal_graded_hubbard8/discovery/boundary_unitary.py',
        'results/marginal_graded_hubbard8/boundary_unitary/certificate.json']
    result={'accepted':True,'seconds':time.monotonic()-start,'rotation':rotation,
            'rotation_family_bound':family,'transferred_targets':targets,
            'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
            'scope':'Exact limit of the proposed one-parameter two-site unitary and physical upper-state transfer to specified nearest-neighbor U,t,V models. Target energies are freshly recomputed from a fixed physical source state. No target ground-energy lower bound or generic molecular transfer is established.'}
    (OUT/'independent_replay.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
