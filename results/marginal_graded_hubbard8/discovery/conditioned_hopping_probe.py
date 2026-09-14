"""Exact spectator-charge-conditioned hopping consistency after spin closure."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_local_hubbard_block import _reflection,_sector
from experiments.marginal_charged_projectors import _particlehole,_spinflip

LABELS=tuple((i,j,k) for i in range(5) for j in range(i+1,5) if (j-i)%2 for k in range(5) if k not in (i,j) and (i,j,k)<(4-j,4-i,4-k))


def operator(label):
    i,j,k=label;weights={}
    for triple,a in [((i,j,k),1),((4-j,4-i,4-k),-1),((i+1,j+1,k+1),-1),((5-j,5-i,5-k),1)]:weights[triple]=weights.get(triple,0)+a
    actions=[]
    for s in range(4096):
        row={}
        for (left,right,spectator),a in weights.items():
            p=((((s>>(2*spectator))&3).bit_count()-1)**2)
            if not a*p:continue
            for spin in (0,1):
                m,n=2*left+spin,2*right+spin
                for word in [((1,m),(0,n)),((1,n),(0,m))]:
                    result=apply_word(word,s)
                    if result:t,sign=result;row[t]=row.get(t,0)+a*p*sign
        actions.append({t:a for t,a in row.items() if a})
    for s,row in enumerate(actions):
        if any(t==s or actions[t].get(s,0)!=a or _sector(s,6)!=_sector(t,6) for t,a in row.items()):raise ValueError('Conditioned hopping physical action fails')
        for transform in [lambda s:_reflection(s,6),_particlehole,_spinflip]:
            us,phase=transform(s);right={}
            for t,a in row.items():
                ut,sign=transform(t);right[ut]=right.get(ut,0)+sign*a
            if right!={t:phase*a for t,a in actions[us].items()}:raise ValueError('Conditioned hopping symmetry fails')
    for sites in (8,9,10):
        totals={}
        for offset in range(sites):
            for (left,right,spectator),a in weights.items():
                pair=tuple(sorted(((left+offset)%sites,(right+offset)%sites)));key=pair+((spectator+offset)%sites,);totals[key]=totals.get(key,0)+a
        if any(totals.values()):raise ValueError('Conditioned hopping translated cancellation fails')
    return actions,2*sum(abs(a) for a in weights.values())


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);folder=parser.parse_args().directory.resolve()
    cp=folder/'range_two_family_limit_certificate.json';rp=folder/'range_two_family_limit_replay.json';r=json.loads(rp.read_text());c=json.loads(cp.read_text())
    if not r['accepted'] or not r.get('spin_telescope') or r['source_sha256'][str(cp.relative_to(ROOT))]!=hashlib.sha256(cp.read_bytes()).hexdigest():raise ValueError('Accepted spin-complete family required')
    results=[]
    if len(LABELS)!=9:raise ValueError('Expected nine canonical odd-distance spectator triples')
    for label in LABELS:
        action,bound=operator(label);moment=F(0)
        for item in c['mixture']:
            v={int(s):a for s,a in item['vector'].items()};norm=sum(a*a for a in v.values());moment+=F(item['weight'])*F(sum(a*b*v.get(t,0) for s,a in v.items() for t,b in action[s].items()),norm)
        results.append({'triple':list(label),'exact_moment':str(moment),'moment_float':float(moment),'norm_upper_bound':bound,'normalized_violation_lower':float(abs(moment)/bound),'violated':bool(moment)})
    files={Path(__file__).resolve(),cp,rp}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'operator':'C(i,j,k)=q_k^2 B(i,j), B=sum_spin(c_i†c_j+h.c.), spectator k distinct from i,j. Y=C(i,j,k)-C(4-j,4-i,4-k); T=Y_left-Y_right.','separators':results,'all_fock_states':4096,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Nine exact necessary translation-consistency tests in a spin-complete physical PSD mixture. Hermiticity, PH/reflection/spin-flip invariance and periodic cancellation checked. Norm upper bound8: spectator q² is a commuting norm-one projector and B has norm2; triangle inequality, not an exact telescope norm. Nonzero moments survive averaging and forbid a stationary quantum extension matching that averaged six-site density matrix. No conditioned-hopping energy certificate or general representability conclusion.'}
    (folder/'conditioned_hopping_overlap.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'accepted':True,'separators':[{k:v for k,v in item.items() if k!='exact_moment'} for item in results]}))


if __name__=='__main__':main()
