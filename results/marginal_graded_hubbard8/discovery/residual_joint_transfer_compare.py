"""Exact same-target transfer of seed, control and enlarged recipes."""
from pathlib import Path
from fractions import Fraction as F
import copy,hashlib,json
from pair_transfer_separation import read_verified
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/marginal_graded_hubbard8/residual_joint'
def read(p):return json.loads(p.read_text())
def main():
    target={'U':'5','t':'1','V':'1/4','W':'-1/5'};files={Path(__file__).resolve(),Path(__file__).with_name('pair_transfer_separation.py')};certs={};receipts={}
    for name in ('seed','control','enlarged'):
        source=BASE/'W_zero'/name/'profile_joint_r1_2_certificate.json';source_r=source.with_name('range_two_replay.json')
        p=BASE/'held_out'/name/'profile_joint_r1_2_certificate.json';rp=p.with_name('range_two_replay.json')
        c=read(source);new=read(p);sr=read_verified(source_r);receipts[name]=read_verified(rp);certs[name]=new;files.update((source,source_r,p,rp))
        if str(source.relative_to(ROOT)) not in sr['source_sha256'] or str(p.relative_to(ROOT)) not in receipts[name]['source_sha256']:raise ValueError('Unbound transfer energy')
        if new['target']!=target:raise ValueError('Wrong held-out target')
        expected=copy.deepcopy(c);old=expected['target'];expected['target']=target;local=expected['local_window']
        local['U']=str(F(5,6)*F(target['U']));local['V']=target['V']
        local['onsite_profile']=[str(F(v)*F(target['U'])/F(old['U'])) for v in local['onsite_profile']]
        local['density_profile']=[str(F(v)*F(target['V'])/F(old['V'])) for v in local['density_profile']]
        local['range_two_density_profile']=[str(F(v)+F(5,4)*(F(target['W'])-F(old['W']))) for v in local['range_two_density_profile']]
        expected['penalized_lower']=new['penalized_lower']
        if expected!=new:raise ValueError('Transfer changed a frozen field')
    p=BASE/'held_out/without_residual/profile_joint_r1_2_certificate.json';rp=p.with_name('range_two_replay.json');c=read(p);receipts['without_residual']=read_verified(rp);files.update((p,rp))
    expected=copy.deepcopy(certs['enlarged']);expected.pop('residual_coherence');expected['kind']='hubbard_projector_extension_v19';expected['penalized_lower']=c['penalized_lower']
    if expected!=c:raise ValueError('Ablation changed old fields')
    if len({r['upper_per_site'] for r in receipts.values()})!=1 or any(r['target']!=target for r in receipts.values()):raise ValueError('Physical upper/target mismatch')
    lo={k:F(r['lower_per_site']) for k,r in receipts.items()}
    result={'accepted':True,'target':target,'open_lower_per_site':{k:str(v) for k,v in lo.items()},
        'enlarged_minus_seed':str(lo['enlarged']-lo['seed']),
        'enlarged_minus_control':str(lo['enlarged']-lo['control']),
        'residual_ablation_contribution':str(lo['enlarged']-lo['without_residual']),
        'physical_upper_per_site':receipts['enlarged']['upper_per_site'],
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
        'scope':'All projector vectors, auxiliary coefficients and penalties frozen from the matched W0 seed, control and enlarged recipes. Only recorded physical profiles rescaled/shifted and scalar threshold recomputed. This one same-geometry/same-filling target does not establish generic transfer or a reoptimized-family ordering.'}
    (BASE/'held_out/comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'accepted':True,**{k:float(F(result[k])) for k in ['enlarged_minus_seed','enlarged_minus_control','residual_ablation_contribution']}}))
if __name__=='__main__':main()
