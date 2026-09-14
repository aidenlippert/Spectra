"""Discover small local witnesses; replay scalable chain intervals exactly."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math,sys,time
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_local_energy_chain import replay_many
OUT=ROOT/'results/marginal_graded_hubbard8/local_energy_chain'


def propose(L,U,upper):
    import numpy as np
    from scipy.linalg import eigh
    data=sector_matrices(L,U,1);lowest=float('inf');best=None;best_energy=float('inf')
    for key,(_,H,columns) in data.items():
        gram=np.diag([sum(a*a for a in column.values()) for column in columns])
        eigenvalues,eigenvectors=eigh(np.array(H,dtype=float),gram,subset_by_index=[0,0])
        e=float(eigenvalues[0]);lowest=min(lowest,e)
        if upper and key[:2]==(L//2,L//2) and e<best_energy:
            best_energy=e;vector={}
            for a,column in zip(eigenvectors[:,0],columns):
                for state,b in column.items():vector[state]=vector.get(state,0)+a*b
            scale=max(abs(a) for a in vector.values())
            best={str(s):int(round(a/scale*10**8)) for s,a in vector.items() if round(a/scale*10**8)}
    certificate={'kind':'local_hubbard_block_v1','sites':L,'U':str(U),'t':'1',
                 'lower':str(F(math.floor(lowest*1000),1000))}
    if upper:certificate['upper_vector']=best
    return certificate,{'sites':L,'U':str(U),'numerical_minimum':lowest,
                        'proposed_lower':certificate['lower'],'all_fock_dimension':4**L,
                        'reflection_blocks':len(data),'maximum_matrix_dimension':max(len(x[1]) for x in data.values())}


def main(replay=False):
    OUT.mkdir(exist_ok=True)
    if replay:
        payload=json.loads((OUT/'certificate.json').read_text())
        chain_certificates=[]
        for entry in payload['chains']:
            chain_certificates.append({'kind':'local_energy_chain_v1','sites':entry['sites'],
                'block_length':entry['block_length'],'U':'4','t':'1',
                'blocks':[payload['local'][name] for name in entry['blocks']],
                'overlap':payload['local'][entry['overlap']]})
        start=time.monotonic();result=replay_many(chain_certificates);result['seconds']=time.monotonic()-start
        result['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in payload['source_sha256']}
        result['source_matches_discovery']=result['source_sha256']==payload['source_sha256']
        (OUT/'independent_replay.json').write_text(json.dumps(result,indent=2)+'\n')
        for row in result['chains']:
            print(json.dumps({'sites':row['sites'],'L':row['block_length'],
                'lower_density':float(F(row['lower_per_site'])),
                'upper_density':float(F(row['upper_per_site'])),
                'width_density':float(F(row['width_per_site']))}),flush=True)
        print(json.dumps({'accepted':True,'seconds':result['seconds'],
                          'unique_local_certificates':result['unique_local_certificates']}),flush=True)
        return
    local={};diagnostics=[]
    for L in (2,4,6):
        for name,U,upper in [(f'physical{L}',F(4),True),(f'overlap{L}',F(4*(L-1),L),False)]:
            certificate,diagnostic=propose(L,U,upper)
            local[name]=certificate;diagnostics.append(diagnostic);print(json.dumps(diagnostic),flush=True)
    chains=[]
    for L in (2,4,6):
        for N in (12,64,1000000):
            remainder=N%L
            chains.append({'sites':N,'block_length':L,'blocks':[f'physical{L}']+([f'physical{remainder}'] if remainder else []),'overlap':f'overlap{L}'})
    sources=['experiments/marginal_local_hubbard_block.py','experiments/marginal_local_energy_chain.py',
             'results/marginal_graded_hubbard8/discovery/local_energy_chain.py',
             'experiments/marginal_polynomial_sos.py','experiments/marginal_transfer_verify.py']
    payload={'local':local,'chains':chains,'discovery':diagnostics,
             'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
             'scope':'Numerical local endpoint/vector proposals only until exact replay. All local sectors must pass. Analytic global composition described in research/marginal_local_energy_chain.md; no whole-chain enumeration.'}
    (OUT/'certificate.json').write_text(json.dumps(payload,indent=2)+'\n')


if __name__=='__main__':main('--replay' in sys.argv)
