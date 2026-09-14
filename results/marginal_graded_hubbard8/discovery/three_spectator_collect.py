"""Collect exact three-spectator energy, family and frozen-transfer evidence."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import re
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/marginal_graded_hubbard8/three_spectator'
OLD=BASE.parent/'two_spectator'


def read(p):return json.loads(p.read_text())
def rel(p):return str(p.relative_to(ROOT))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n')


def main():
    log=(BASE/'full_validation.log').read_text()
    count=re.search(r'(\d+) passed',log);subtests=re.search(r'(\d+) subtests passed',log)
    if not count or ' failed' in log:raise ValueError('Completed passing full suite required')
    suites=list(ET.parse(BASE/'full_validation.xml').getroot().iter('testsuite'))
    if not suites or any(int(s.attrib.get(k,0)) for s in suites for k in ('failures','errors')):
        raise ValueError('Passing JUnit required')
    if '58 passed' not in (BASE/'focused_validation.log').read_text():raise ValueError('Focused production gates missing')
    if '3 passed' not in (BASE/'fraction_free_validation.log').read_text():raise ValueError('Fraction-free tests missing')
    paths=[];source_files=set();entries=0;construction_entries=0
    for p in BASE.rglob('*.json'):
        if p in (BASE/'summary.json',BASE/'provenance.json'):continue
        d=read(p)
        if not isinstance(d,dict):continue
        for source,h in d.get('source_sha256',{}).items():
            if sha(ROOT/source)!=h:raise ValueError('Stale current construction/proof source: '+source)
            source_files.add(ROOT/source);construction_entries+=1
        if d.get('accepted') is True and 'source_sha256' in d:
            paths.append(p);entries+=len(d['source_sha256'])
    cases={}
    names=('range_two_replay','range_two_family_limit_replay','previous_family_comparison',
           'three_spectator_closure','two_spectator_closure','one_spectator_overlap',
           'pair_transfer_overlap','spin_overlap','coherent_overlap','charge_markov_extension')
    for case in ('W_zero','W_plus_1'):
        d=BASE/case/'final'
        for name in names:
            if read(d/(name+'.json')).get('accepted') is not True:raise ValueError('Missing selected proof '+case+'/'+name)
        energy=read(d/'range_two_replay.json');family=read(d/'range_two_family_limit_replay.json')
        comparison=read(d/'previous_family_comparison.json')
        if not family.get('three_spectator_hopping'):raise ValueError('Enlarged family mode missing')
        three=read(d/'three_spectator_closure.json')['three_spectator_moments']
        two=read(d/'two_spectator_closure.json')['two_spectator_moments']
        if len(three)!=18 or any(map(F,three.values())) or len(two)!=30 or any(map(F,two.values())):
            raise ValueError('Exact enlarged closures missing')
        for name in ('pair_transfer_overlap','spin_overlap'):
            if any(F(item['exact_moment']) for item in read(d/(name+'.json'))['separators']):raise ValueError('Old closure failed')
        if any(map(F,read(d/'one_spectator_overlap.json')['spectator_moments'].values())) or F(read(d/'coherent_overlap.json')['exact_moment']):
            raise ValueError('Old hopping closure failed')
        lower=F(energy['lower_replay']['periodic_lower_density']);cap=F(family['periodic_family_upper']);gap=cap-lower
        if gap<0 or gap!=F(family['family_gap']) or lower!=F(family['accepted_periodic_lower']):raise ValueError('Exact family interval mismatch')
        cases[case]={'directory':rel(d),'periodic_lower':str(lower),'family_ceiling':str(cap),
                     'family_gap':str(gap),'family_gap_float':float(gap),
                     'open_lower_per_site':energy['lower_per_site'],'physical_upper_per_site':energy['upper_per_site'],
                     'mixture_sources':family['family_replay']['nearest_constraint_replay']['mixture_sources'],
                     'lower_improvement':comparison['lower_improvement'],
                     'strict_previous_family_separation':comparison['strict_family_separation'],
                     'signed_previous_family_separation':comparison['exact_separation'],
                     'signed_previous_family_separation_float':comparison['separation_float'],
                     'three_spectator_zero_moments':18,'two_spectator_zero_moments':30}
    # Preserve the distinction between current replays and earlier source versions.
    prior=read(OLD/'summary.json');refined=read(OLD/'refinement_summary.json')
    historical=refined['current_receipts']+prior['current_receipts']+prior['historical_audit']['receipt_paths']
    if len(historical)!=len(set(historical)):raise ValueError('Duplicate historical receipts')
    snapshot_dirs=[BASE.parent/'pair_transfer/source_before',BASE.parent/'pair_transfer/family_source_before',OLD/'source_before',BASE/'source_before']
    snapshots={}
    for directory in snapshot_dirs:
        for name,h in read(directory/'original_sha256.json').items():
            p=directory/Path(name).name
            if sha(p)!=h:raise ValueError('Altered original source snapshot')
            snapshots[(name,h)]=p
    unchanged=archived=0
    for name in historical:
        d=read(ROOT/name)
        if d.get('accepted') is not True:raise ValueError('Unaccepted historical receipt')
        for source,h in d['source_sha256'].items():
            if sha(ROOT/source)==h:unchanged+=1
            elif (source,h) in snapshots:archived+=1
            else:raise ValueError('Unaccounted historical source change '+source)
    files=source_files|set(paths)|{Path(__file__).resolve(),ROOT/'research/marginal_three_spectator_energy.md',OLD/'summary.json',OLD/'refinement_summary.json'}
    files.update(p for p in BASE.rglob('*') if p.is_file() and p not in (BASE/'summary.json',BASE/'provenance.json'))
    for directory in snapshot_dirs:files.update(p for p in directory.rglob('*') if p.is_file())
    files.update(ROOT/name for name in historical)
    discovery=BASE.parent/'discovery'
    files.update(discovery.glob('three_spectator*.py'))
    for name in ('joint_congruence_replay.py','signed_charge_numeric.py','joint_profile_numeric.py',
                 'two_spectator_prune_ledger.py','two_spectator_family_trust.py','two_spectator_family_resume.py',
                 'two_spectator_primal_enrich.py','pair_transfer_separation.py','two_spectator_overlap_replay.py'):
        files.add(discovery/name)
    files.update(ROOT/'tests'/name for name in ('test_marginal_three_spectator_hopping.py',
                 'test_marginal_three_spectator_matching.py','test_marginal_three_spectator_fraction_free.py'))
    diagnostics={}
    for p in BASE.rglob('*.json'):
        if p.name not in ('diagonal_family_limit_diagnostic.json','thermal_proposal.json','subspace_diagnostic.json','weight_limit_refusal.json'):continue
        d=read(p)
        diagnostics[rel(p)]={k:v for k,v in d.items() if k in ('accepted','proposal_written','matrix_evaluations',
            'optimizer_success','optimizer_message','seconds','states','mixture_sources','rational_family_upper_float',
            'proposed_periodic_lower_float','basis_size','exact_basis_failure','refusal_observed','error',
            'maximum_weight_characters','new_vector_integer_scale','numerical_objective','status','history')}
        if 'pricing' in d:diagnostics[rel(p)].update(pricing_rounds=len(d['pricing']),last_pricing=d['pricing'][-1:])
    transfer=read(BASE/'held_out/frozen_transfer_comparison.json')
    if transfer.get('accepted') is not True:raise ValueError('Frozen comparison missing')
    summary={'accepted':True,'cases':cases,'frozen_transfer':transfer,'current_receipts':list(map(rel,paths)),
             'current_receipts_audited':len(paths),'current_source_hash_entries_checked':entries,
             'construction_source_hash_entries_checked':construction_entries,
             'historical_audit':{'receipts':len(historical),'unchanged_hash_entries':unchanged,
                 'preserved_snapshot_entries':archived,'receipt_paths':historical,
                 'scope':'Historical receipt sources audited against preserved snapshots where current source changed; not fresh current-verifier replay.'},
             'full_validation':{'tests_passed':int(count[1]),'subtests_passed':int(subtests[1]),
                 'junit_suites':[s.attrib for s in suites],
                 'command':'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q --junitxml=results/marginal_graded_hubbard8/three_spectator/full_validation.xml',
                 'scope':'Final production changes and matching tests included. Three later fraction-free discovery tests passed separately.'},
             'focused_tests_passed':58,'additional_fraction_free_tests_passed':3,
             'energy_version':16,'family_version':12,'source_cap':137,'diagnostics':diagnostics,
             'report':'research/marginal_three_spectator_energy.md','provenance':rel(BASE/'provenance.json'),
             'remaining':['Both enlarged family gaps remain nonzero; no exact numerical-limit or attainment claim.',
                          'A nonpositive signed comparison does not establish separation from the preceding entire family.',
                          'Tested stationary moments vanish; full quantum extension has not been established.',
                          'Frozen coupling transfer only. General representability, generic molecular/long-range/higher-dimensional transfer and requested-accuracy scalability remain unproved. Goal active.']}
    write(BASE/'provenance.json',{'scope':'Current exact proof and construction provenance. Numeric proposals and refused candidates remain nonaccepting. Historical source snapshots retained.',
                                  'sha256':{rel(p):sha(p) for p in sorted(files)}})
    write(BASE/'summary.json',summary)
    central=ROOT/'results/marginal_final_validation.json';d=read(central);d['three_spectator_energy']=summary
    d['latest_validation_scope']='ENERGYv16/FAMILYv12 integrate18 three-spectator conditions with exact enlarged moments and unchanged old-mode gates. Stronger exact lower certificates, enlarged family caps and frozen transfer/ablation accepted. Signed preceding-family comparisons distinguish proved from unproved separation. Full production regression plus3 later discovery tests pass. Numerical limits and general quantum representability remain open.'
    write(central,d)
    print(json.dumps({'accepted':True,'current_receipts':len(paths),'current_source_entries':entries,
                      'historical_receipts':len(historical),'historical_unchanged_entries':unchanged,
                      'historical_snapshot_entries':archived,'construction_entries':construction_entries,
                      'provenance_files':len(files),'full_tests':int(count[1]),'subtests':int(subtests[1])}))


if __name__=='__main__':main()
