"""Exact spin-dot translation consistency after signed-charge/hopping closure."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,sys
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_local_hubbard_block import _reflection,_sector
from experiments.marginal_charged_projectors import _particlehole,_spinflip


def spin_pair_terms(i,j):
    # 4 S_i.S_j = (n_up-n_down)_i(n_up-n_down)_j
    #             + 2(S_i^+ S_j^- + S_i^- S_j^+).
    terms=[]
    for a in (0,1):
        for b in (0,1):
            m,n=2*i+a,2*j+b
            terms.append(((((1,m),(0,m),(1,n),(0,n))),(-1)**(a+b)))
    for a in (0,1):
        terms.append(((((1,2*i+a),(0,2*i+1-a),(1,2*j+1-a),(0,2*j+a))),2))
    return terms


def operator(i,j):
    weights={}
    for pair,a in [((i,j),1),((4-j,4-i),-1),((i+1,j+1),-1),((5-j,5-i),1)]:weights[pair]=weights.get(pair,0)+a
    words=[(word,a*b) for (left,right),a in weights.items() if a for word,b in spin_pair_terms(left,right)]
    action=[]
    for s in range(4096):
        image={}
        for word,a in words:
            result=apply_word(word,s)
            if result:
                t,sign=result;image[t]=image.get(t,0)+a*sign
        action.append({t:a for t,a in image.items() if a})
    for s,image in enumerate(action):
        for t,a in image.items():
            if action[t].get(s,0)!=a or _sector(t,6)!=_sector(s,6):raise ValueError('Spin operator fails Hermiticity or spin-number conservation')
        for transform in [lambda s:_reflection(s,6),_particlehole,_spinflip]:
            us,phase=transform(s);other={}
            for t,a in image.items():
                ut,sign=transform(t);other[ut]=other.get(ut,0)+sign*a
            if other!={t:phase*a for t,a in action[us].items()}:raise ValueError('Spin operator symmetry failure')
    for sites in (8,9,10):
        totals={}
        for offset in range(sites):
            for (left,right),a in weights.items():
                pair=tuple(sorted(((left+offset)%sites,(right+offset)%sites)));totals[pair]=totals.get(pair,0)+a
        if any(totals.values()):raise ValueError('Translated spin telescope fails')
    # Each pair operator is zero outside two singly occupied sites and has
    # singlet/triplet eigenvalues -3,+1. Confirm its full two-site CAR action.
    pair=[[0]*16 for _ in range(16)]
    for s in range(16):
        for word,a in spin_pair_terms(0,1):
            result=apply_word(word,s)
            if result:t,sign=result;pair[t][s]+=a*sign
    expected=[[0]*16 for _ in range(16)]
    expected[5][5]=expected[10][10]=1
    expected[6][6]=expected[9][9]=-1
    expected[6][9]=expected[9][6]=2
    if pair!=expected:raise ValueError('Independent two-spin singlet/triplet block fails')
    return action,3*sum(abs(a) for a in weights.values())


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);folder=parser.parse_args().directory.resolve()
    cp=folder/'range_two_family_limit_certificate.json';rp=folder/'range_two_family_limit_replay.json';r=json.loads(rp.read_text())
    if not r['accepted'] or not r.get('hopping_telescope') or r['source_sha256'][str(cp.relative_to(ROOT))]!=hashlib.sha256(cp.read_bytes()).hexdigest():raise ValueError('Accepted hopping-complete family required')
    c=json.loads(cp.read_text());results=[]
    for i,j in [(0,1),(0,2),(0,3),(1,2)]:
        action,bound=operator(i,j);moment=F(0)
        for item in c['mixture']:
            v={int(s):a for s,a in item['vector'].items()};norm=sum(a*a for a in v.values());moment+=F(item['weight'])*F(sum(a*b*v.get(t,0) for s,a in v.items() for t,b in action[s].items()),norm)
        results.append({'pair':[i,j],'exact_moment':str(moment),'moment_float':float(moment),'norm_upper_bound':bound,'normalized_violation_lower':float(abs(moment)/bound),'violated':bool(moment)})
    files={Path(__file__).resolve(),cp,rp}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'operator_convention':'A(i,j)=4 S_i dot S_j; Y=A(i,j)-A(4-j,4-i); T=Y_left-Y_right','separators':results,'all_fock_states':4096,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Exact necessary translation-consistency constraints for a stationary quantum extension. PH/reflection/spin-flip invariant CAR actions, all-Fock Hermiticity and periodic cancellation verified. Pair norm3 from exact singlet/triplet block; triangle upper bound12 for each telescope, not a claim of its exact norm. No spin-corrected energy bound or general representability result.'}
    (folder/'spin_overlap.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({**result,'separators':[{k:v for k,v in item.items() if k!='exact_moment'} for item in results],'source_sha256':'omitted'}))


if __name__=='__main__':main()
