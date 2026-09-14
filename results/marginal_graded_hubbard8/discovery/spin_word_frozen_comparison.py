"""Exact frozen transfer, two-term ablation and previous-recipe comparison."""
from pathlib import Path
from fractions import Fraction as F
import copy,hashlib,json
from pair_transfer_separation import read_verified

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/marginal_graded_hubbard8/spin_word'

def main():
    source=BASE/'W_zero/profile_joint_r1_2_certificate.json'
    paths=[BASE/'held_out'/name/'profile_joint_r1_2_certificate.json' for name in ('','without_coherent','previous')]
    original=json.loads(source.read_text());new,ablation,previous=[json.loads(p.read_text()) for p in paths]
    source_receipt=source.with_name('range_two_replay.json');read_verified(source_receipt)
    receipts=[p.with_name('range_two_replay.json') for p in paths]
    energy,no_new,old=[read_verified(p) for p in receipts]
    expected=copy.deepcopy(original);target=new['target'];old_target=original['target']
    if target!={'U':'5','t':'1','V':'1/4','W':'-1/5'}:raise ValueError('Explicit held-out target required')
    expected['target']=target;local=expected['local_window']
    local['U']=str(F(5,6)*F(target['U']))
    local['onsite_profile']=[str(F(v)*F(target['U'])/F(old_target['U'])) for v in local['onsite_profile']]
    local['V']=target['V'];local['density_profile']=[str(F(v)*F(target['V'])/F(old_target['V'])) for v in local['density_profile']]
    local['range_two_density_profile']=[str(F(v)+F(5,4)*(F(target['W'])-F(old_target['W']))) for v in local['range_two_density_profile']]
    expected['penalized_lower']=new['penalized_lower']
    if expected!=new:raise ValueError('Frozen transfer changed a correction or physical recipe')
    expected=copy.deepcopy(new);expected.pop('spin_word_telescope');expected['kind']='hubbard_projector_extension_v17';expected['penalized_lower']=ablation['penalized_lower']
    if expected!=ablation:raise ValueError('Ablation changed other coefficients')
    previous_source=BASE.parent/'coherent_projector/held_out/profile_joint_r1_2_certificate.json'
    if previous!=json.loads(previous_source.read_text()):raise ValueError('Previous frozen recipe changed')
    if len({r['upper_per_site'] for r in (energy,no_new,old)})!=1 or any(r['target']!=target for r in (energy,no_new,old)):
        raise ValueError('Physical upper or target mismatch')
    lower=F(energy['lower_per_site']);loss=lower-F(no_new['lower_per_site']);improvement=lower-F(old['lower_per_site'])
    files=[Path(__file__).resolve(),Path(__file__).with_name('pair_transfer_separation.py'),source,source_receipt,previous_source]+paths+receipts
    result={'accepted':True,'target':target,'source_recipe':'Initial W0 spin-word certificate, before polish.',
            'with_spin_word_open_lower':str(lower),'without_spin_word_open_lower':no_new['lower_per_site'],
            'previous_frozen_open_lower':old['lower_per_site'],'exact_ablation_loss':str(loss),'ablation_loss_float':float(loss),
            'exact_previous_recipe_improvement':str(improvement),'previous_recipe_improvement_float':float(improvement),
            'physical_upper':energy['upper_per_site'],
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            'scope':'Exact frozen-source/penalty/correction transfer to U5,t1,V1/4,W-1/5 on the same geometry and filling. Full diagonal fixed-recipe ablation and previous frozen-recipe comparison use accepted energy proofs. Neither comparison caps a reoptimized older family. No generic molecular or scalability conclusion.'}
    (BASE/'held_out/frozen_transfer_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'accepted':True,'ablation_loss':float(loss),'previous_recipe_improvement':float(improvement)}))

if __name__=='__main__':main()
