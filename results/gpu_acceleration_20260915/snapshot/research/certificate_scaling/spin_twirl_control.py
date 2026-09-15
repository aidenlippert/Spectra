"""Source-certificate SU(2) feasibility control, explicitly NOT H-only discovery."""
from fractions import Fraction as F
from pathlib import Path
from collections import defaultdict
import argparse,hashlib,json,time
import numpy as np
from scipy import sparse
from experiments.marginal_coefficient import dictionaries
from experiments.marginal_symbolic import canonical,mono,decode,encode,verify
from research.certificate_scaling.spin_basis import decompose_words,rref
from research.certificate_scaling.spin_parity import spin_parity
from research.certificate_scaling.spin_twirl import twirl
from research.certificate_scaling.spin_invariant_discovery import export


def coordinate_maps(groups):
    """Invert each constant-size local pattern exactly; map monomials to copies."""
    patterns=defaultdict(list)
    for gi,g in enumerate(groups):
        for a,copy in enumerate(g['copies']):
            for r,p in enumerate(copy):
                w=next(iter(p));key=(tuple(i//2 for c,i in w if c),tuple(i//2 for c,i in w if not c))
                patterns[key].append(((gi,a,r),p))
    coordinates={};max_local=0
    for entries in patterns.values():
        words=sorted({w for _,p in entries for w in p});n=len(words);max_local=max(max_local,n)
        if n!=len(entries):raise ValueError('Local spin basis must be square')
        A=[[p.get(w,F(0)) for _,p in entries] for w in words]
        inv,pivots=rref([row+[F(i==j) for j in range(n)] for i,row in enumerate(A)],2*n)
        if pivots!=list(range(n)):raise ValueError('Dependent local spin basis')
        inverse=[row[n:] for row in inv]
        for j,w in enumerate(words):
            coordinates[w]=[(key,inverse[i][j]) for i,(key,_) in enumerate(entries) if inverse[i][j]]
            # Explicit coefficient reconstruction, independently of pivot bookkeeping.
            for k,u in enumerate(words):
                if sum(A[k][i]*inverse[i][j] for i in range(n))!=F(k==j):raise AssertionError('Bad exact coordinate inverse')
    return coordinates,max_local


def run(source,out):
    start=time.monotonic();raw=Path(source).read_bytes();cert=json.loads(raw);m=cert['modes'];n=cert['particles'];out=Path(out);out.mkdir(parents=True,exist_ok=False)
    h=decode(cert['hamiltonian'],m,4);masks=spin_parity(h,m);groups=[];families=[]
    for family in dictionaries(m,'mixed'):
        words=[w for w in family['words'] if len(w)!=1]
        if not words:continue
        gs,_=decompose_words(words,m,masks)
        for g in gs:g['family']=family['name']
        offset=len(groups);mapping,maxlocal=coordinate_maps(gs)
        mapping={w:[((gi+offset,a,r),v) for (gi,a,r),v in terms] for w,terms in mapping.items()}
        groups.extend(gs);families.append((family['name'],mapping,maxlocal))
    grams=[np.zeros((len(g['copies']),len(g['copies']))) for g in groups]
    block_matches=[];term_updates=0
    for block in cert['blocks']:
        normalized=[canonical(mono(tuple(tuple(x) for x in w))) for w in block['words']]
        if any(len(p)!=1 for p in normalized):raise ValueError('Source dictionary is not monomial')
        words=[next(iter(p)) for p in normalized]
        matches=[(name,mapping,local) for name,mapping,local in families if set(words)<=mapping.keys()]
        if len(matches)!=1:raise ValueError('Source block is outside the unambiguous degree-2/3 spin family')
        name,mapping,_=matches[0];block_matches.append(name)
        Y=np.asarray(block['factor'],dtype=float)/cert['denominator']
        if Y.ndim!=2 or Y.shape[1]!=len(words) or not np.isfinite(Y).all():raise ValueError('Invalid source factors')
        transforms=defaultdict(lambda:([],[],[]))
        for col,(w,p) in enumerate(zip(words,normalized)):
            sign=p[w]
            for (gi,a,r),value in mapping[w]:
                rows,cols,values=transforms[gi,r];rows.append(col);cols.append(a);values.append(float(value*sign));term_updates+=1
        for (gi,r),(rows,cols,values) in transforms.items():
            g=groups[gi];T=sparse.csr_matrix((values,(rows,cols)),shape=(len(words),len(g['copies'])))
            C=Y@T;grams[gi]+=(C.T@C)/((g['two_spin']+1)*g['weights'][r])
    transformed=time.monotonic()
    np.savez_compressed(out/'twirled_grams.npz',**{f'q{i}':q for i,q in enumerate(grams)})
    X=twirl(decode(cert['number_multiplier'],m,4))
    X={w:c for w,c in X.items() if all(sum((mask>>i)&1 for _,i in w)%2==0 for mask in masks)}
    blocks,den,negative=export(groups,grams)
    result_cert={**cert,'hamiltonian':encode(h),'number_multiplier':encode(X),'blocks':blocks,'denominator':den}
    exported=time.monotonic();exact=verify(result_cert);data=json.dumps(result_cert,separators=(',',':'))+'\n';(out/'certificate.json').write_text(data)
    result={'source_certificate':str(source),'source_sha256':hashlib.sha256(raw).hexdigest(),'source_certificate_used':True,
        'scope':'Feasibility control transformed from an existing SOS. This does not discover a certificate from H alone; original source discovery and independent upper costs remain.',
        'modes':m,'particles':n,'source_blocks':len(cert['blocks']),'block_families':block_matches,'Gram_entries':sum(q.size for q in grams),
        'coordinate_nonzero_updates':term_updates,'max_local_coordinate_dimension':max(x[2] for x in families),
        'transform_seconds':transformed-start,'export_seconds':exported-transformed,'exact_replay_seconds':time.monotonic()-exported,
        'wall_seconds':time.monotonic()-start,'negative_eigenvalue_mass':negative,'certificate_bytes':len(data.encode()),'exact':exact}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='block_families'}),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.source,a.out)
