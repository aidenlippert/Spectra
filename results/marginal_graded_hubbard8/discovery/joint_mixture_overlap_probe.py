"""Discover and exactly check a violated five-site overlap-consistency entry."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import _reflection,_sector
from experiments.marginal_charged_projectors import _transform,_spinflip,_particlehole,charged_vectors
OUT=ROOT/'results/marginal_graded_hubbard8/joint_projector/signed_density'


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--symmetry-average',action='store_true');args=parser.parse_args()
    source=OUT/'family_limit_certificate.json';c=json.loads(source.read_text())
    states=[(F(item['weight']),{int(s):a for s,a in item['vector'].items()}) for item in c['mixture']]
    if args.symmetry_average:
        half={int(s):a for s,a in c['half_vector'].items() if a};charged,_=charged_vectors(c['charged_vector'])
        for fn in (_spinflip,_particlehole):
            if not any(_transform(half,fn)=={s:p*a for s,a in half.items()} for p in (-1,1)):
                raise ValueError('Half projector not symmetry invariant')
            for v in charged:
                if not any(_transform(v,fn)=={s:p*a for s,a in u.items()} for p in (-1,1) for u in charged):
                    raise ValueError('Charged family not symmetry invariant')
        states=[(w/4,u) for w,v in states for u in (v,_transform(v,_spinflip),_transform(v,_particlehole),_transform(_transform(v,_spinflip),_particlehole))]
    delta=np.zeros((1024,1024))
    for w,v in states:
        norm=sum(a*a for a in v.values())
        left=np.zeros((4,1024));right=np.zeros((4,1024))
        for s,a in v.items():left[s>>10,s&1023]=a;right[s&3,s>>2]=a
        delta+=float(w/norm)*(left.T@left-right.T@right)
    a,b=np.unravel_index(np.argmax(abs(delta)),delta.shape);a,b=int(a),int(b)
    def entry(i,j):
        total=F(0)
        for w,v in states:
            norm=sum(x*x for x in v.values())
            total+=w*F(sum(v.get(i+1024*z,0)*v.get(j+1024*z,0)-v.get(4*i+z,0)*v.get(4*j+z,0) for z in range(4)),norm)
        return total
    exact=entry(a,b)
    if not exact or _sector(a,5)!=_sector(b,5):raise ValueError('No conserved nonzero overlap separator')
    if args.symmetry_average:
        def spin5(s):return sum((((s>>(2*i))&1)*2+((s>>(2*i+1))&1))<<(2*i) for i in range(5))
        visited=set();candidates=[]
        for s in range(1024):
            if s in visited:continue
            orbit={s,spin5(s),1023^s,1023^spin5(s)}
            reflected={_reflection(t,5)[0] for t in orbit};visited.update(orbit|reflected)
            if orbit&reflected:continue
            Ydiag={**{t:1 for t in orbit},**{t:-1 for t in reflected}}
            expectation=sum(v*entry(t,t) for t,v in Ydiag.items())
            candidates.append({'diagonal':{str(t):v for t,v in sorted(Ydiag.items())},'exact_expectation':str(expectation),'expectation_float':float(expectation)})
        candidates.sort(key=lambda c:abs(c['expectation_float']),reverse=True)
        (OUT/'symmetry_diagonal_candidates.json').write_text(json.dumps({'basis_size':len(candidates),'candidates':candidates,'scope':'Reflection-odd, spin-exchange and particle-hole invariant five-site diagonal orbit basis; exact moments in symmetry-averaged family-limit mixture.'},indent=2)+'\n')
    # Antisymmetrize the Hermitian matrix unit under five-site reflection.
    Y={(a,b):1}
    if a!=b:Y[(b,a)]=1
    for (i,j),v in list(Y.items()):
        ri,si=_reflection(i,5);rj,sj=_reflection(j,5)
        Y[(ri,rj)]=Y.get((ri,rj),0)-si*sj*v
    Y={k:v for k,v in Y.items() if v}
    expectation=sum(v*entry(j,i) for (i,j),v in Y.items())
    reflected={}
    for (i,j),v in Y.items():
        ri,si=_reflection(i,5);rj,sj=_reflection(j,5);reflected[(ri,rj)]=si*sj*v
    if reflected!={k:-v for k,v in Y.items()} or not expectation:raise ValueError('Invalid reflection-odd separator')
    result={'accepted':True,'mixture_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'symmetry_averaged':args.symmetry_average,
            'selected_entry':[a,b],'exact_overlap_difference':str(exact),'difference_float':float(exact),
            'numerical_frobenius_norm':float(np.linalg.norm(delta)),
            'maximum_diagonal_difference':float(np.max(abs(np.diag(delta)))),
            'reflection_odd_five_site_matrix':[[i,j,v] for (i,j),v in sorted(Y.items())],
            'exact_telescoping_expectation':str(expectation),'telescoping_expectation_float':float(expectation),
            'scope':'Exact nonzero five-site marginal-overlap separator for the local family-limit mixture. Such a mixture is not a translation-invariant extendible local state. The reflection-odd five-site Y gives a reflection-even six-site telescoping correction Y_left-Y_right; its periodic translated sum is zero. No improved energy certificate is claimed here.'}
    (OUT/('symmetry_averaged_overlap_separator.json' if args.symmetry_average else 'overlap_consistency_separator.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('exact_overlap_difference','exact_telescoping_expectation')}),flush=True)


if __name__=='__main__':main()
