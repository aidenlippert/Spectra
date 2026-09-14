"""Exact PH-even one-spectator charge-conditioned hopping basis and rank."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_local_hubbard_block import _reflection
from experiments.marginal_charged_projectors import _particlehole,_spinflip,charged_vectors
from experiments.marginal_polynomial_sos import integer_psd

LABELS=tuple((i,j,k,2 if (j-i)%2 else 1) for i in range(5) for j in range(i+1,5) for k in range(5) if k not in (i,j) and (i,j,k)<(4-j,4-i,4-k))


def operator(label):
    i,j,k,power=label;weights={}
    for triple,a in [((i,j,k),1),((4-j,4-i,4-k),-1),((i+1,j+1,k+1),-1),((5-j,5-i,5-k),1)]:weights[triple]=weights.get(triple,0)+a
    action=[]
    for s in range(4096):
        row={}
        for (left,right,spectator),a in weights.items():
            value=((((s>>(2*spectator))&3).bit_count()-1)**power)
            for spin in (0,1):
                m,n=2*left+spin,2*right+spin
                for word in [((1,m),(0,n)),((1,n),(0,m))]:
                    result=apply_word(word,s)
                    if result:t,sign=result;row[t]=row.get(t,0)+a*value*sign
        action.append({t:a for t,a in row.items() if a})
    for s,row in enumerate(action):
        if any((s^t).bit_count()!=2 or action[t].get(s,0)!=a for t,a in row.items()):raise ValueError('Physical single-hop Hermiticity failed')
        for transform in [lambda s:_reflection(s,6),_particlehole,_spinflip]:
            us,phase=transform(s);right={}
            for t,a in row.items():
                ut,sign=transform(t);right[ut]=right.get(ut,0)+sign*a
            if right!={t:phase*a for t,a in action[us].items()}:raise ValueError('Physical symmetry failed')
    for sites in (8,9,10):
        totals={}
        for offset in range(sites):
            for (left,right,spectator),a in weights.items():
                pair=tuple(sorted(((left+offset)%sites,(right+offset)%sites)));key=pair+((spectator+offset)%sites,) if power else pair;totals[key]=totals.get(key,0)+a
        if any(totals.values()):raise ValueError('Translated telescope failed')
    return action


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);args=parser.parse_args();folder=args.directory.resolve()
    cp=folder/'range_two_family_limit_certificate.json';rp=folder/'range_two_family_limit_replay.json';r=json.loads(rp.read_text());c=json.loads(cp.read_text())
    if not r['accepted'] or not r.get('spin_telescope') or r['source_sha256'][str(cp.relative_to(ROOT))]!=hashlib.sha256(cp.read_bytes()).hexdigest():raise ValueError('Accepted spin family required')
    if len(LABELS)!=14 or sum(label[-1]==2 for label in LABELS)!=9:raise ValueError('Expected nine quadratic and five linear spectator classes')
    charged,_=charged_vectors(c['charged_vector']);particles={int(s).bit_count() for s,a in c['half_vector'].items() if a}|{s.bit_count() for vector in charged for s,a in vector.items() if a}
    if particles!={5,6,7}:raise ValueError('Previous fixed source particle sectors differ')
    old=[operator(label) for label in [(0,1,2,0),(1,2,0,0),(0,3,1,0)]];new=[operator(label) for label in LABELS];data=old+new
    gram=[[sum(a*right[s].get(t,0) for s,row in enumerate(left) if s.bit_count()<=2 for t,a in row.items()) for right in data] for left in data];rank=integer_psd(gram)
    if rank['rank']!=17:raise ValueError('Conditioned span has unexpected dependencies')
    results=[]
    oldprobe=json.loads((folder/'conditioned_hopping_overlap.json').read_text())
    for label,action in zip(LABELS,new):
        moment=F(0)
        for item in c['mixture']:
            v={int(s):a for s,a in item['vector'].items()};norm=sum(a*a for a in v.values());moment+=F(item['weight'])*F(sum(a*b*v.get(t,0) for s,a in v.items() for t,b in action[s].items()),norm)
        if label[-1]==2:
            previous=next(item for item in oldprobe['separators'] if item['triple']==list(label[:3]))
            if F(previous['exact_moment'])!=moment:raise ValueError('Nine existing quadratic probes disagree')
        results.append({'triple':list(label[:3]),'charge_power':label[-1],'exact_moment':str(moment),'moment_float':float(moment),'norm_upper_bound':8,'violated':bool(moment)})
    files={Path(__file__).resolve(),cp,rp,folder/'conditioned_hopping_overlap.json'}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'separators':results,'gram':gram,'exact_rank':rank,'source_particle_numbers':sorted(particles),'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Complete nonconstant one-spectator charge-function hopping class after PH-even/reflection-odd projection: q^2 for odd-distance hopping (nine), q for even-distance hopping (five). Unconditioned odd-distance hopping is the three old directions. Rank17 of old3+new14 on N<=2 proves independence modulo old diagonal and spin corrections (zero single-hop entries), and fixed projectors (N5/6/7). Does not claim completeness of all quantum consistency operators. Exact CAR symmetry, periodic cancellation and moments checked; conservative norm bound8.'}
    out=ROOT/'results/marginal_graded_hubbard8/spectator_hopping'/folder.parent.name;out.mkdir(parents=True,exist_ok=True);(out/'basis_probe.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'accepted':True,'rank':rank['rank'],'separators':[{k:v for k,v in item.items() if k!='exact_moment'} for item in results]}))


if __name__=='__main__':main()
