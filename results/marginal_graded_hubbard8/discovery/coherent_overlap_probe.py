"""Exact off-diagonal translation-consistency separator after charge closure."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,sys
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_local_hubbard_block import _reflection
from experiments.marginal_charged_projectors import _particlehole,_spinflip


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);folder=parser.parse_args().directory.resolve()
    cp=folder/'range_two_family_limit_certificate.json';rp=folder/'range_two_family_limit_replay.json';receipt=json.loads(rp.read_text())
    if not receipt['accepted'] or not receipt.get('full_signed_charge') or receipt['source_sha256'][str(cp.relative_to(ROOT))]!=hashlib.sha256(cp.read_bytes()).hexdigest():raise ValueError('Accepted complete signed-charge family required')
    words=[]
    for i,weight in enumerate([1,-2,1]):
        for spin in range(2):
            a,b=2*i+spin,2*(i+3)+spin
            words.extend([(((1,a),(0,b)),weight),(((1,b),(0,a)),weight)])
    action=[]
    for s in range(4096):
        row={}
        for word,weight in words:
            result=apply_word(word,s)
            if result:
                target,sign=result;row[target]=row.get(target,0)+weight*sign
        action.append({target:a for target,a in row.items() if a})
    for s,row in enumerate(action):
        if row.get(s,0):raise ValueError('Expected a purely off-diagonal operator')
        if any(action[target].get(s,0)!=a for target,a in row.items()):raise ValueError('Hermiticity fails')
        for transform in [lambda state:_reflection(state,6),_particlehole,_spinflip]:
            us,phase=transform(s);left={target:phase*a for target,a in action[us].items()};right={}
            for target,a in row.items():
                ut,sign=transform(target);right[ut]=right.get(ut,0)+sign*a
            if left!=right:raise ValueError('Physical symmetry invariance fails')
    # The six disjoint mode-pair hoppings have norm one each, so the weighted
    # triangle bound is eight. Construct an exact +8 eigenvector to attain it.
    eigenvector={0:1}
    for spin in range(2):
        for i,sign in enumerate([1,-1,1]):
            new={}
            for s,a in eigenvector.items():
                for mode,coefficient in [(2*i+spin,1),(2*(i+3)+spin,sign)]:
                    target,phase=apply_word(((1,mode),),s);new[target]=new.get(target,0)+a*coefficient*phase
            eigenvector={s:a for s,a in new.items() if a}
    image={}
    for s,a in eigenvector.items():
        for target,b in action[s].items():image[target]=image.get(target,0)+a*b
    if {s:a for s,a in image.items() if a}!={s:8*a for s,a in eigenvector.items()}:raise ValueError('Norm-attaining eigenvector fails')
    for sites in (8,9,10):
        totals={}
        for offset in range(sites):
            for i,a in enumerate([1,-2,1]):
                pair=tuple(sorted(((offset+i)%sites,(offset+i+3)%sites)));totals[pair]=totals.get(pair,0)+a
        if any(totals.values()):raise ValueError('Periodic translated sum fails')
    c=json.loads(cp.read_text());moment=F(0)
    for item in c['mixture']:
        vector={int(s):a for s,a in item['vector'].items()};norm=sum(a*a for a in vector.values())
        numerator=sum(a*b*vector.get(target,0) for s,a in vector.items() for target,b in action[s].items());moment+=F(item['weight'])*F(numerator,norm)
    files={Path(__file__).resolve(),cp,rp}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'operator':'B(0,3)-2B(1,4)+B(2,5), B(i,j)=sum_spin(c_i^dagger*c_j+c_j^dagger*c_i)','five_site_telescope':'Y=B(0,3)-B(1,4)','exact_moment':str(moment),'moment_float':float(moment),'local_operator_norm':8,'normalized_violation':float(abs(moment)/8),'violated':bool(moment),'car_words':len(words),'full_fock_symmetry_checks':4096,'norm_attaining_vector':{str(s):a for s,a in eigenvector.items()},'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Exact noncommuting translation-consistency separator in an accepted charge-complete local PSD mixture. Hermiticity, PH/reflection/spin-exchange invariance and norm8 checked; its translated sum vanishes. A nonzero moment persists under PH/reflection averaging and forbids a stationary quantum extension matching that averaged local density matrix. No stronger energy certificate or general representability claim is made.'}
    (folder/'coherent_overlap.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','exact_moment','norm_attaining_vector')}))


if __name__=='__main__':main()
