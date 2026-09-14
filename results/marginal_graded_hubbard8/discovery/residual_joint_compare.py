"""Source-bound matched continuation controls and exact frozen-recipe ceilings."""
from pathlib import Path
from fractions import Fraction as F
import copy,hashlib,json,sys
from pair_transfer_separation import read_verified
from residual_coherence_fixed_limit_replay import replay as fixed_replay
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/marginal_graded_hubbard8/residual_joint'
FIELDS=('penalty','penalized_lower','spin_word_telescope','hopping_telescope','spin_telescope','spectator_hopping','pair_transfer','two_spectator_hopping','three_spectator_hopping','coherent_projector','pure_coherence','residual_coherence','kind')
def read(p):return json.loads(p.read_text())
def fixed(c):
    c=copy.deepcopy(c)
    for k in FIELDS:c.pop(k,None)
    c['local_window'].pop('hopping_profile');c['joint'].pop('penalty')
    return c

def hashes(files):
    files.add(Path(__file__).resolve())
    for module in tuple(sys.modules.values()):
        n=getattr(module,'__file__',None)
        if n:
            p=Path(n).resolve()
            if p.is_relative_to(ROOT/'experiments') or p.is_relative_to(BASE.parent/'discovery'):files.add(p)
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}

def main():
    for case,selection in [('W_zero','refined'),('W_plus_1','scaled')]:
        files=set();certs={};energies={};families={}
        for name in ('seed','control','enlarged','without_residual'):
            p=BASE/case/name/'profile_joint_r1_2_certificate.json';rp=p.with_name('range_two_replay.json')
            certs[name]=read(p);energies[name]=read_verified(rp);files.update((p,rp))
            if str(p.relative_to(ROOT)) not in energies[name]['source_sha256']:raise ValueError('Energy certificate binding missing')
            if name in ('seed','control'):
                fp=p.with_name('range_two_family_limit_replay.json');families[name]=read_verified(fp);files.add(fp)
        seed,control,new,no_new=[certs[k] for k in ('seed','control','enlarged','without_residual')]
        if seed['kind']!='hubbard_projector_extension_v19' or control['kind']!=seed['kind'] or new['kind']!='hubbard_projector_extension_v20':raise ValueError('Explicit matched versions required')
        if not fixed(seed)==fixed(control)==fixed(new):raise ValueError('Joint search changed fixed family data')
        expected=copy.deepcopy(new);expected.pop('residual_coherence');expected['kind']=seed['kind'];expected['penalized_lower']=no_new['penalized_lower']
        if no_new!=expected:raise ValueError('Ablation changed old coefficients')
        original=BASE.parent/'pure_coherence'/case/selection/'profile_joint_r1_2_certificate.json'
        if original.read_bytes()!=(BASE/case/'seed/profile_joint_r1_2_certificate.json').read_bytes():raise ValueError('Seed recipe changed')
        proposal_path=BASE.parent/'residual_coherence'/case/'fixed_limit/family_proposal.json';proposal=read(proposal_path)
        if (ROOT/proposal['seed_certificate']).read_bytes()!=original.read_bytes():raise ValueError('Frozen coefficient ceiling uses another seed')
        cap=fixed_replay(seed,proposal);files.update((original,proposal_path))
        lower={k:F(v['lower_replay']['periodic_lower_density']) for k,v in energies.items()}
        if F(cap['accepted_seed_lower'])!=lower['seed']:raise ValueError('Fixed ceiling seed lower mismatch')
        family_ceiling=F(families['seed']['periodic_family_upper'])
        if F(families['control']['periodic_family_upper'])!=family_ceiling:raise ValueError('Control does not share old family ceiling')
        if len({v['upper_per_site'] for v in energies.values()})!=1:raise ValueError('Physical upper mismatch')
        result={'accepted':True,'case':case,'periodic_lower':{k:str(v) for k,v in lower.items()},
            'gain_over_seed':str(lower['enlarged']-lower['seed']),
            'gain_over_equal_budget_control':str(lower['enlarged']-lower['control']),
            'gain_over_best_old_certificate':str(lower['enlarged']-max(lower['seed'],lower['control'])),
            'fixed_recipe_ablation_contribution':str(lower['enlarged']-lower['without_residual']),
            'preceding_family_ceiling':str(family_ceiling),
            'signed_preceding_full_family_separation':str(lower['enlarged']-family_ceiling),
            'strict_preceding_full_family_separation':lower['enlarged']>family_ceiling,
            'fresh_fixed_coefficient_cap':cap,
            'joint_excess_above_frozen_ceiling':str(lower['enlarged']-F(cap['fixed_recipe_family_ceiling'])),
            'source_sha256':hashes(files),
            'scope':'Exact energy comparison of matched bounded continuations and a fixed-recipe ablation. Fresh physical ceiling caps only the two new coefficients with every seed field frozen. Neither numerical continuation is an optimizer-attainment proof; matched-control superiority does not establish whole-family separation.'}
        out=BASE/case/'comparison.json';out.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({'case':case,'accepted':True,**{k:float(F(result[k])) for k in ['gain_over_seed','gain_over_equal_budget_control','joint_excess_above_frozen_ceiling','signed_preceding_full_family_separation']}}))
if __name__=='__main__':main()
