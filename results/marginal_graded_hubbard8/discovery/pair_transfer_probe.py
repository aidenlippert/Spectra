"""Exact pair-transfer consistency after one-spectator hopping closure."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_local_hubbard_block import _reflection,_sector
from experiments.marginal_charged_projectors import _particlehole,_spinflip,charged_vectors
from experiments.marginal_polynomial_sos import integer_psd

PAIRS=((0,1),(0,2),(0,3),(1,2))


def pair_terms(i,j):
    return [((1,2*i),(1,2*i+1),(0,2*j+1),(0,2*j)),((1,2*j),(1,2*j+1),(0,2*i+1),(0,2*i))]


def operator(i,j):
    weights={}
    for pair,a in [((i,j),1),((4-j,4-i),-1),((i+1,j+1),-1),((5-j,5-i),1)]:weights[pair]=weights.get(pair,0)+a
    action=[]
    for s in range(4096):
        row={}
        for (left,right),a in weights.items():
            for word in pair_terms(left,right):
                result=apply_word(word,s)
                if result:t,sign=result;row[t]=row.get(t,0)+a*sign
        action.append({t:a for t,a in row.items() if a})
    for s,row in enumerate(action):
        if any((s^t).bit_count()!=4 or _sector(s,6)!=_sector(t,6) or action[t].get(s,0)!=a for t,a in row.items()):raise ValueError('Pair-transfer action fails')
        for transform in [lambda s:_reflection(s,6),_particlehole,_spinflip]:
            us,phase=transform(s);right={}
            for t,a in row.items():
                ut,sign=transform(t);right[ut]=right.get(ut,0)+sign*a
            if right!={t:phase*a for t,a in action[us].items()}:raise ValueError('Pair-transfer symmetry fails')
    for sites in (8,9,10):
        totals={}
        for offset in range(sites):
            for (left,right),a in weights.items():
                pair=tuple(sorted(((left+offset)%sites,(right+offset)%sites)));totals[pair]=totals.get(pair,0)+a
        if any(totals.values()):raise ValueError('Pair-transfer periodic cancellation fails')
    return action


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);folder=parser.parse_args().directory.resolve();cp=folder/'range_two_family_limit_certificate.json';rp=folder/'range_two_family_limit_replay.json';r=json.loads(rp.read_text());c=json.loads(cp.read_text())
    if not r['accepted'] or not r.get('spectator_hopping') or r['source_sha256'][str(cp.relative_to(ROOT))]!=hashlib.sha256(cp.read_bytes()).hexdigest():raise ValueError('Accepted full one-spectator family required')
    for s in range(16):
        row={}
        for word in pair_terms(0,1):
            result=apply_word(word,s)
            if result:t,sign=result;row[t]=row.get(t,0)+sign
        expected={12:1} if s==3 else ({3:1} if s==12 else {})
        if row!=expected:raise ValueError('Independent two-site pair transfer block fails')
    data=[operator(*pair) for pair in PAIRS];doubles=[3<<(2*i) for i in range(6)]
    gram=[[sum(a*right[s].get(t,0) for s in doubles for t,a in left[s].items()) for right in data] for left in data];rank=integer_psd(gram)
    if rank['rank']!=4:raise ValueError('Four independent pair directions required')
    charged,_=charged_vectors(c['charged_vector']);particles={int(s).bit_count() for s,a in c['half_vector'].items() if a}|{s.bit_count() for v in charged for s,a in v.items() if a}
    if particles!={5,6,7}:raise ValueError('Previous fixed projector sectors differ')
    results=[]
    for pair,action in zip(PAIRS,data):
        moment=F(0)
        for item in c['mixture']:
            v={int(s):a for s,a in item['vector'].items()};norm=sum(a*a for a in v.values());moment+=F(item['weight'])*F(sum(a*b*v.get(t,0) for s,a in v.items() for t,b in action[s].items()),norm)
        results.append({'pair':list(pair),'exact_moment':str(moment),'moment_float':float(moment),'norm_upper_bound':4,'normalized_violation_lower':float(abs(moment)/4),'violated':bool(moment)})
    files={Path(__file__).resolve(),cp,rp}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'operator':'J(i,j)=d_i† d_j+d_j† d_i, d_i†=c_i_up† c_i_down†. Y=J(i,j)-J(4-j,4-i), T=Y_left-Y_right.','separators':results,'double_occupancy_gram':gram,'exact_rank':rank,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Four exact stationary quantum consistency tests. Full-Fock CAR Hermiticity and PH/reflection/spin-flip symmetry, periodic cancellation and independent norm-one two-site pair-transfer block checked; telescope triangle norm bound4. Rank4 on the two-electron double-occupancy subspace proves independence modulo the previous family: diagonal terms have no off-diagonal entries, spin exchange cannot move a double, single hopping changes two rather than four bits, and fixed projectors occupy N5/6/7. Nonzero moments survive averaging and forbid a stationary extension of that averaged six-site matrix. No pair-corrected energy or general representability claim.'}
    (folder/'pair_transfer_overlap.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'accepted':True,'rank':4,'separators':[{k:v for k,v in item.items() if k!='exact_moment'} for item in results]}))


if __name__=='__main__':main()
