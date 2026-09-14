"""Read-only exact certificate verification, independent of the exporter.

CLI emits a receipt to stdout. Certificate data are never rewritten.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import copy, hashlib, json, sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_determinant_tree import DeterminantOracle

HERE = Path(__file__).resolve().parent
EXPECTED = '78365286838362d1b591e7afcef0f1e2a085fafc445e961f3b7e19603dbf1c6c'
def require(condition, message):
    if not condition: raise ValueError(message)

def dot(a,b): return sum(x*b.get(s,0) for s,x in a.items())
def gram(a,b): return [[dot(x,y) for y in b] for x in a]
def matchings(sites):
    if not sites: return [()]
    result=[]
    for k in range(1,len(sites),2):
        for inside in matchings(sites[1:k]):
            for outside in matchings(sites[k+1:]):
                result.append(((sites[0],sites[k]),)+inside+outside)
    return result

def rank(a):
    a=[[F(x) for x in row] for row in a]; r=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][j]),None)
        if pivot is None: continue
        a[r],a[pivot]=a[pivot],a[r]; scale=a[r][j]
        a[r]=[x/scale for x in a[r]]
        for i in range(r+1,len(a)):
            scale=a[i][j]
            a[i]=[x-scale*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a): break
    return r

def ldl(a):
    n=len(a)
    require(all(len(row)==n for row in a),'square matrix')
    require(all(a[i][j]==a[j][i] for i in range(n) for j in range(n)),'symmetric matrix')
    L=[[F(i==j) for j in range(n)] for i in range(n)]; d=[]
    for j in range(n):
        d.append(F(a[j][j])-sum(L[j][k]**2*d[k] for k in range(j)))
        require(d[-1]!=0,'nonzero unpivoted LDL pivot')
        for i in range(j+1,n):
            L[i][j]=(a[i][j]-sum(L[i][k]*L[j][k]*d[k] for k in range(j)))/d[j]
    require(all(sum(L[i][k]*d[k]*L[j][k] for k in range(n))==a[i][j]
                for i in range(n) for j in range(n)), 'LDL exact reconstruction')
    return d

def verify(c):
    h=json.loads((HERE.parent/'hamiltonian.json').read_text())
    digest=hashlib.sha256(json.dumps(h,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    require(digest==EXPECTED==c['hamiltonian_sha256'],'original Hamiltonian digest')
    require(h['modes']==16 and h['particles']==8,'original sector')
    require(c['schema']=='h8-singlet-valence-embedding-v1' and c['sites']==8,'schema/site count')
    p=[sum(1<<(2*i+(i not in up)) for i in range(8)) for up in combinations(range(8),4)]
    require(c['valence_states']==p,'complete valence coordinates')
    ms=matchings(tuple(range(8)))
    require(c['matchings']==[[list(pair) for pair in m] for m in ms] and c['matching_count']==14,'complete noncrossing matchings')
    # Reconstruct each coefficient from its spin assignment and site permutation.
    expected=[]
    for pairs in ms:
        order=[i for pair in pairs for i in pair]
        phase=(-1)**sum(order[i]>order[j] for i in range(8) for j in range(i+1,8))
        v={}
        for s in p:
            spins=[(s>>(2*i))&1 for i in range(8)]
            if all(spins[i]!=spins[j] for i,j in pairs):
                v[s]=phase*(-1)**sum(not spins[i] for i,j in pairs)
        expected.append(v)
    require(type(c['basis']) is list and len(c['basis'])==14,'basis count')
    V=[]
    for data in c['basis']:
        require(type(data) is dict and all(type(a) is int for a in data.values()),'integer basis')
        V.append({int(s):a for s,a in data.items()})
    require(V==expected,'basis coefficient reconstruction')
    targets=set(); plus=[]
    for s in p:
        v={}
        for i in range(8):
            z=apply_word(((1,2*i),(0,2*i+1)),s)
            if z: v[z[0]]=v.get(z[0],0)+z[1]
        plus.append(v);targets.update(v)
    for v in V:
        require(all(sum(v.get(s,0)*plus[k].get(t,0) for k,s in enumerate(p))==0 for t in targets),'S+ annihilation')
    actual_rank=rank([[col.get(t,0) for col in plus] for t in sorted(targets)])
    require(actual_rank==56==c['splus_rank'] and len(targets)==56==c['splus_target_dimension'],'S+ exact rank')
    require(70-actual_rank==14==c['splus_kernel_dimension'],'complete singlet kernel')
    G=gram(V,V); gd=ldl(G)
    require(G==c['gram'] and all(x>0 for x in gd) and list(map(str,gd))==c['gram_ldl'],'positive exact Gram')
    oracle=DeterminantOracle(h); actions=[oracle.action(s) for s in p]
    require(all(row.get(s,0)==0 for row in actions for s in p),'actual H_PP zero')
    W=gram(actions,actions)
    require([[F(x) for x in row] for row in c['leakage_gram']]==W,'actual leakage Gram')
    images=[]
    for v in V:
        out={}
        for k,s in enumerate(p):
            for t,a in actions[k].items():out[t]=out.get(t,0)+v.get(s,0)*a
        images.append({s:a for s,a in out.items() if a})
    WG=gram(images,images)
    require([[F(x) for x in row] for row in c['projected_leakage']]==WG,'projected leakage')
    gamma,tau=F(c['gamma']),F(c['tau'])
    require(gamma==F(-381,100) and tau==F(-21,5),'fixed diagnostic thresholds')
    S=[[-tau*G[i][j]-WG[i][j]/(gamma-tau) for j in range(14)] for i in range(14)]
    sd=ldl(S); negative=sum(x<0 for x in sd)
    require(list(map(str,sd))==c['restricted_schur_ldl'] and negative==13==c['restricted_schur_negative_pivots'],'exact Schur inertia')
    return {'checks':'passed','hamiltonian_sha256':digest,'splus_annihilation':True,'gram_rank':14,
            'splus_rank':actual_rank,'actual_h_pp_zero':True,'negative_pivots':negative,
            'unique_action_states':len(oracle.cache),'referenced_determinants':oracle.referenced_state_count(),
            'scope':'Exact complete valence singlet embedding and fixed scalar-Schur inertia only; no ground energy certificate or verification of Lieb theorem.'}

def main():
    c=json.loads((HERE/'certificate.json').read_text()); receipt=verify(c)
    changes={
        'duplicate_basis':lambda q:q['basis'].__setitem__(1,copy.deepcopy(q['basis'][0])),
        'altered_coefficient':lambda q:q['basis'][0].__setitem__(next(iter(q['basis'][0])),99),
        'wrong_gram':lambda q:q['gram'][0].__setitem__(0,999),
        'wrong_spin_rank':lambda q:q.__setitem__('splus_rank',55),
        'wrong_schur_pivot':lambda q:q['restricted_schur_ldl'].__setitem__(0,'1'),
        'wrong_threshold':lambda q:q.__setitem__('gamma','-3'),
    }
    results={}
    for name,change in changes.items():
        q=copy.deepcopy(c);change(q)
        try: verify(q)
        except ValueError as error: results[name]=str(error)
        else: raise ValueError('Tampered certificate was accepted: '+name)
    receipt['tamper_rejections']=results
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
