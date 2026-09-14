"""Numerical profile discovery followed by exact family and global replay."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_window_family_bound import replay as replay_family
from experiments.marginal_local_energy_chain import replay_many
OUT=ROOT/'results/marginal_graded_hubbard8/weighted_window_family'
BASE=ROOT/'results/marginal_graded_hubbard8/local_energy_chain/certificate.json'


def discover():
    import numpy as np
    from scipy.linalg import eigh
    from scipy.optimize import minimize
    def data(a,b):return sector_matrices(4,3,1,[a,6-a,6-a,a],[b,3-2*b,b])
    ds=[data(0,0),data(1,0),data(0,1)];family=[]
    for key,(_,H,cols) in ds[0].items():
        base=np.array(H,float);A=np.array(ds[1][key][1],float)-base;B=np.array(ds[2][key][1],float)-base
        G=np.diag([sum(a*a for a in c.values()) for c in cols])
        family.append((key,base,A,B,G,cols))
    calls=0
    def evaluate(x):
        nonlocal calls
        calls+=1;best=None
        for key,H,A,B,G,cols in family:
            values,vectors=eigh(H+x[0]*A+x[1]*B,G,subset_by_index=[0,0]);v=vectors[:,0]
            if best is None or values[0]<best[0]:best=(float(values[0]),np.array([v@A@v,v@B@v]),key,v,cols)
        return best
    def objective(x):
        value,gradient,*_=evaluate(x)
        return -value,-gradient
    proposal=minimize(objective,[.53,.75],jac=True,bounds=[(0,6),(0,1.5)],method='L-BFGS-B',options={'ftol':1e-15,'gtol':1e-11,'maxiter':100,'maxfun':150})
    value,gradient,key,v,cols=evaluate(proposal.x);physical={}
    for amplitude,column in zip(v,cols):
        for state,sign in column.items():physical[state]=physical.get(state,0)+amplitude*sign
    scale=max(abs(a) for a in physical.values())
    vector={str(s):int(round(a/scale*10**12)) for s,a in physical.items() if round(a/scale*10**12)}
    certificate={'kind':'hubbard_four_window_family_bound_v1','a':'531373/1000000','b':'3/4',
                 'lower':'-2.040424675','upper_vector':vector,'target_family_upper':'-2.040424674'}
    diagnostic={'numerical_parameters':proposal.x.tolist(),'numerical_minimum':value,
                'numerical_gradient':gradient.tolist(),'minimizing_sector':list(key),'evaluations':calls,
                'optimizer_message':str(proposal.message),'scope':'Numerical proposal only; exact replay supplies acceptance and the family-wide bound.'}
    return certificate,diagnostic


def main(replay=False):
    OUT.mkdir(exist_ok=True)
    if not replay:
        certificate,diagnostic=discover()
        (OUT/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
        (OUT/'discovery.json').write_text(json.dumps(diagnostic,indent=2)+'\n')
        print(json.dumps(diagnostic),flush=True)
        return
    certificate=json.loads((OUT/'certificate.json').read_text());start=time.monotonic()
    family=replay_family(certificate)
    (OUT/'family_independent_replay.json').write_text(json.dumps(family,indent=2)+'\n')
    print(json.dumps({k:family[k] for k in ('local_minimum_lower','maximum_local_minimum_upper','family_bracket_width','affine_expectation')}),flush=True)
    old=json.loads(BASE.read_text());chains=[]
    # Replay the original nine intervals under the current source as well.
    for entry in old['chains']:
        chains.append({'kind':'local_energy_chain_v1','sites':entry['sites'],'block_length':entry['block_length'],
                       'U':'4','t':'1','blocks':[old['local'][name] for name in entry['blocks']],
                       'overlap':old['local'][entry['overlap']]})
    for N in (12,64,1000000):
        remainder=N%6
        blocks=[old['local']['physical6']]+([old['local'][f'physical{remainder}']] if remainder else [])
        chains.append({'kind':'local_energy_chain_v1','sites':N,'block_length':6,'U':'4','t':'1',
                       'blocks':blocks,'overlap':family['local_certificate']})
    result=replay_many(chains);result['seconds_including_family']=time.monotonic()-start
    sources=['experiments/marginal_local_hubbard_block.py','experiments/marginal_local_energy_chain.py',
             'experiments/marginal_window_family_bound.py','results/marginal_graded_hubbard8/discovery/weighted_window_family.py']
    result['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}
    (OUT/'global_independent_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in result['chains'][-3:]:
        print(json.dumps({'sites':row['sites'],'lower_density':float(F(row['lower_per_site'])),
                          'upper_density':float(F(row['upper_per_site'])),'width_density':float(F(row['width_per_site']))}),flush=True)
    print(json.dumps({'accepted':True,'seconds':result['seconds_including_family'],'unique_local_certificates':result['unique_local_certificates']}),flush=True)

if __name__=='__main__':main('--replay' in sys.argv)
