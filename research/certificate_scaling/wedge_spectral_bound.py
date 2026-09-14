"""Rational Gram witnesses for exterior-power residual eigenvalue bounds.

NumPy proposes factors; the independent replay uses integer arithmetic only.
The indices are k-tuples, never individual orbitals for k>1.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, lcm
from pathlib import Path
import argparse, hashlib, json, sys, time
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import verified_residual
from research.certificate_scaling.wedge_residual_bound import residual_bounds

MAX_COMPONENT = 600
MAX_CUBIC_WORK = 150_000_000


def matrix_lower(matrix, witness):
    n=len(matrix)
    if not n or any(len(row)!=n for row in matrix): raise ValueError('Square nonempty matrix required')
    matrix=[[F(v) for v in row] for row in matrix]
    if any(matrix[i][j]!=matrix[j][i] for i in range(n) for j in range(i)): raise ValueError('Hermitian matrix required')
    factor=witness['factor']; denominator=witness['denominator']; ell=F(witness['ell'])
    if type(denominator) is not int or denominator<=0 or type(factor) is not list or len(factor)!=n:
        raise ValueError('Integer factors and positive denominator required')
    rank=len(factor[0])
    if any(type(row) is not list or len(row)!=rank or any(type(x) is not int for x in row) for row in factor):
        raise ValueError('Invalid integer factor shape')
    if n>MAX_COMPONENT or rank>n or n*n*rank>MAX_CUBIC_WORK: raise ValueError('Wedge proof cost gate')
    den=lcm(denominator**2,ell.denominator,*(v.denominator for row in matrix for v in row))
    gram_scale=den//(denominator**2); shift=int(ell*den)
    diagonal=[int(matrix[i][i]*den)-shift-sum(x*x for x in factor[i])*gram_scale for i in range(n)]
    row_abs=[0]*n
    for i in range(n):
        for j in range(i):
            residual=int(matrix[i][j]*den)-sum(x*y for x,y in zip(factor[i],factor[j]))*gram_scale
            row_abs[i]+=abs(residual); row_abs[j]+=abs(residual)
    correction=F(min(diagonal[i]-row_abs[i] for i in range(n)),den)
    return ell+correction, {'dimension':n,'rank':rank,'integer_product_work':n*(n+1)//2*rank,
        'ell':str(ell),'residual_gershgorin':str(correction),'lower':str(ell+correction)}


def body_components(residual,m,k):
    indices=list(combinations(range(m),k)); adjacency={i:set() for i in indices}; entries={}
    phase=(-1)**(k*(k-1)//2)
    for word,c in residual.items():
        if len(word)!=2*k:continue
        if [cr for cr,_ in word]!=[1]*k+[0]*k:raise ValueError('Balanced canonical body required')
        i=tuple(p for _,p in word[:k]);j=tuple(p for _,p in word[k:])
        if i not in adjacency or j not in adjacency:raise ValueError('Noncanonical wedge index')
        entries[i,j]=c*phase
        if i!=j and c:adjacency[i].add(j);adjacency[j].add(i)
    if any(entries.get((j,i),F(0))!=c for (i,j),c in entries.items()):raise ValueError('Non-Hermitian body')
    unseen=set(indices);components=[]
    while unseen:
        first=min(unseen);unseen.remove(first);stack=[first];component=[]
        while stack:
            i=stack.pop();component.append(i)
            for j in adjacency[i]:
                if j in unseen:unseen.remove(j);stack.append(j)
        component.sort();components.append(component)
    return [(part,[[entries.get((i,j),F(0)) for j in part] for i in part]) for part in components]


def extract(cert):
    degree=cert.get('operator_degree',3)
    if degree!=3:raise ValueError('Cubic residual scope required')
    return verified_residual(cert)


def propose(residual,m,n):
    import numpy as np
    bodies=[];work=0
    for k in range(1,min(3,n)+1):
        parts=[]
        for indices,matrix in body_components(residual,m,k):
            d=len(indices);work+=d**3
            if d>MAX_COMPONENT or work>MAX_CUBIC_WORK:raise ValueError('Total spectral proposal cost gate')
            ev,u=np.linalg.eigh(np.array(matrix,dtype=float))
            ell=F(int(np.floor((float(ev[0])-1e-12)*10**15)),10**15)
            numeric=u*np.sqrt(np.maximum(ev-float(ell),0))[None,:]
            if not np.all(np.isfinite(numeric)) or np.max(np.abs(numeric),initial=0)>1e6:raise ValueError('Factor numeric range')
            factor=np.rint(numeric*10**12).astype(np.int64).tolist()
            parts.append({'indices':[list(i) for i in indices],'ell':str(ell),'factor':factor,'denominator':10**12})
        bodies.append({'body':k,'components':parts})
    return {'method':'fixed_N_wedge_spectral_gram_v1','bodies':bodies,'proposal_cubic_work':work}


def replay_residual(cert,residual,witness):
    m,n=cert['modes'],cert['particles']; old_correction,old_details=residual_bounds(residual,m,n)
    expected=list(range(1,min(3,n)+1))
    if [b['body'] for b in witness['bodies']]!=expected:raise ValueError('Missing/repeated body proof')
    lower=F(cert['b'])+residual.get((),F(0));details=[];work=0
    for body in witness['bodies']:
        k=body['body'];parts=body_components(residual,m,k)
        if len(parts)!=len(body['components']):raise ValueError('Missing component proof')
        endpoints=[];receipts=[]
        for (indices,matrix),proof in zip(parts,body['components']):
            if proof['indices']!=[list(i) for i in indices]:raise ValueError('Wrong exterior-power indices')
            work+=len(indices)**3
            if work>MAX_CUBIC_WORK:raise ValueError('Total exact proof cost gate')
            endpoint,receipt=matrix_lower(matrix,proof);endpoints.append(endpoint);receipts.append(receipt)
        spectral=min(endpoints,default=F(0))*comb(n,k)
        previous=next((F(b['chosen_lower']) for b in old_details if b['body']==k),F(0))
        chosen=max(previous,spectral);lower+=chosen
        details.append({'body':k,'wedge_dimension':comb(m,k),'components':receipts,'spectral_fixed_N_lower':str(spectral),'chosen_lower':str(chosen)})
    if lower<F(cert['b'])+old_correction:raise AssertionError('Weaker than existing residual bound')
    return {'method':witness['method'],'lower':str(lower),'lower_float':float(lower),'bodies':details,
        'old_wedge_lower':str(F(cert['b'])+old_correction),'many_body_states_enumerated':0,
        'scope':'Exact residual bound only; full original SOS discovery cost remains.'}


def run(path,out,proof_path=None):
    start=time.monotonic();raw=path.read_bytes();cert=json.loads(raw);residual,old=extract(cert)
    extracted=time.monotonic()
    if proof_path:
        witness=json.loads(proof_path.read_text())
        if witness['certificate_sha256']!=hashlib.sha256(raw).hexdigest():raise ValueError('Certificate hash mismatch')
    else:
        witness=propose(residual,cert['modes'],cert['particles']);witness['certificate_sha256']=hashlib.sha256(raw).hexdigest()
    proposed=time.monotonic();receipt=replay_residual(cert,residual,witness)
    receipt.update({'coefficient_lower':old['lower'],'certificate_sha256':hashlib.sha256(raw).hexdigest(),
        'extract_seconds':extracted-start,'proposal_seconds':proposed-extracted,'spectral_replay_seconds':time.monotonic()-proposed,'wall_seconds':time.monotonic()-start})
    out.mkdir(parents=True,exist_ok=True)
    if not proof_path:(out/'witness.json').write_text(json.dumps(witness,separators=(',',':'))+'\n')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='bodies'}),flush=True)
    return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--proof',type=Path);a=p.parse_args();run(a.certificate,a.out,a.proof)
