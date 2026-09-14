"""Exact cyclic-state ceiling for every six-site weighting of one target.

Optional --discover proposes the physical ring state numerically. The default
path uses only standard-library exact arithmetic and validates its witness.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_local_hubbard_block import _terms
OUT=ROOT/'results/marginal_graded_hubbard8/density_transfer'

def rotate(state):
    out={}
    for s,a in state.items():
        moved=[(i+2)%12 for i in range(12) if s>>i&1]
        inversions=sum(x>y for j,x in enumerate(moved) for y in moved[j+1:])
        out[sum(1<<i for i in moved)]=a*(-1)**inversions
    return out

def observables(state):
    q=[];hopping=[];density=[]
    for i in range(6):
        j=(i+1)%6
        q.append(sum(a*a*(F(((s>>(2*i))&3)==3)-F(((s>>(2*i))&3).bit_count()-1,2)) for s,a in state.items()))
        density.append(sum(a*a*(((s>>(2*i))&3).bit_count()-1)*(((s>>(2*j))&3).bit_count()-1) for s,a in state.items()))
        value=0
        for s,a in state.items():
            for spin in (0,1):
                left,right=2*i+spin,2*j+spin
                for word in (((1,left),(0,right)),((1,right),(0,left))):
                    image=apply_word(word,s)
                    if image:value+=a*image[1]*state.get(image[0],0)
        hopping.append(value)
    return q,hopping,density

def five_site_rdm(state,left):
    groups={}
    for s,a in state.items():
        local,environment=(s&1023,s>>10) if left else (s>>2,s&3)
        groups.setdefault(environment,{})[local]=a
    matrix={}
    for vector in groups.values():
        for i,a in vector.items():
            for j,b in vector.items():matrix[i,j]=matrix.get((i,j),0)+a*b
    return {ij:a for ij,a in matrix.items() if a}

def discover():
    import numpy as np
    from scipy.linalg import eigh
    states=[s for s in range(4096) if s.bit_count()==6 and sum((s>>(2*i))&1 for i in range(6))==3]
    ix={s:i for i,s in enumerate(states)};matrix=np.zeros((400,400))
    terms=_terms(6,4,1,density=[F(1,2)]*5)
    for spin in (0,1):
        a,b=10+spin,spin;terms.extend([(((1,a),(0,b)),-1),(((1,b),(0,a)),-1)])
    for j,s in enumerate(states):
        matrix[j,j]+=float(F(1,2)*((s&3).bit_count()-1)*(((s>>10)&3).bit_count()-1))
        for word,a in terms:
            image=apply_word(word,s)
            if image:matrix[ix[image[0]],j]+=float(a*image[1])
    values,vectors=eigh(matrix,subset_by_index=[0,0]);vector=vectors[:,0];vector/=max(abs(vector))
    state={str(s):int(round(float(a)*10**9)) for s,a in zip(states,vector) if int(round(float(a)*10**9))}
    (OUT/'window_ceiling_certificate.json').write_text(json.dumps({'sites':6,'target':{'U':'4','t':'1','V':'1/2'},
        'physical_state':state,'claimed_profile_density_ceiling':'-0.66302383'},indent=2)+'\n')

def main():
    if '--discover' in sys.argv:discover()
    c=json.loads((OUT/'window_ceiling_certificate.json').read_text())
    if c['sites']!=6 or c['target']!={'U':'4','t':'1','V':'1/2'}:raise ValueError('Specified six-site target required')
    raw=c['physical_state']
    if type(raw) is not dict or not 1<=len(raw)<=4096:raise ValueError('Bounded physical state required')
    state={}
    for key,a in raw.items():
        s=int(key)
        if s in state or not 0<=s<4096 or s.bit_count()!=6 or type(a) is not int or abs(a)>10**9:
            raise ValueError('Invalid physical six-particle amplitude')
        state[s]=a
    norm=sum(a*a for a in state.values())
    if norm<=0:raise ValueError('Nonzero physical norm required')
    sums=[[0]*6 for _ in range(3)];current=state;initial=observables(state)
    for _ in range(6):
        if sum(a*a for a in current.values())!=norm:raise ValueError('Translation changed norm')
        for group,values in zip(sums,observables(current)):
            for i,value in enumerate(values):group[i]+=value
        current=rotate(current)
    if current!=state:raise ValueError('Six fermionic translations did not close')
    means=[]
    for group in sums:
        if len(set(group))!=1:raise ValueError('Cyclic mixture has nonuniform local observables')
        means.append(F(group[0],6*norm))
    q,h,d=means
    family_density=4*q-h+F(1,2)*d
    ring=F(4*sum(initial[0])-sum(initial[1])+F(1,2)*sum(initial[2]),norm)
    if family_density!=ring/6:raise ValueError('Independent ring-to-window counting disagrees')
    claimed=F(c['claimed_profile_density_ceiling'])
    if family_density>claimed:raise ValueError('Physical witness does not certify the stated ceiling')
    # The submitted witness is itself exactly translation invariant, so
    # its local consistency obstruction can also be checked without mixing.
    if rotate(state)!=state:raise ValueError('Pure cyclic witness required for this extension obstruction')
    left,right=five_site_rdm(state,True),five_site_rdm(state,False)
    if left!=right:raise ValueError('Five-site overlaps do not agree entrywise')
    purity=F(sum(a*a for a in left.values()),norm*norm)
    if not 0<purity<1:raise ValueError('Mixed five-site reduction required')
    sources=['results/marginal_graded_hubbard8/discovery/density_window_ceiling.py',
        'results/marginal_graded_hubbard8/density_transfer/window_ceiling_certificate.json',
        'experiments/marginal_transfer_verify.py','experiments/marginal_local_hubbard_block.py']
    result={'accepted':True,'target':c['target'],'physical_support':len(state),'norm':str(norm),
        'uniform_centered_onsite':str(q),'uniform_positive_hopping':str(h),'uniform_density_correlation':str(d),
        'ring_trial_energy':str(ring),'profile_lower_ceiling':str(5*family_density),
        'profile_density_ceiling':str(family_density),'outward_profile_density_ceiling':str(claimed),
        'pure_cyclic_state_exact':True,'five_site_marginals_equal_entrywise':True,
        'five_site_purity':str(purity),'repeated_seven_site_extension_excluded':True,
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        'scope':'Ceiling on ALL six-site onsite/hopping/density profile lower bounds with sums20,5,5/2, also allowing arbitrary Hermitian even five-site telescoping boundary corrections. The physical cyclic state makes the profile expectations identical and the boundary-correction expectation zero. Its two five-site marginals agree but are mixed; purity of the full six-site state excludes a seven-site state with that same pure marginal on both consecutive six-site windows. No profile sign or reflection restriction is needed for the ceiling. This finite-ring trial energy is NOT a bulk ground-energy bound.'}
    (OUT/'window_ceiling_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Accepted exact profile density ceiling:',float(family_density),flush=True)

if __name__=='__main__':main()
