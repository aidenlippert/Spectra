"""Audit the refined family cap and new independent three-spectator obstructions."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8/two_spectator'
DIRECTORIES = ('trust_seed','trust','trust_refined','trust_final','primal_subspace',
               'primal_atoms','primal_completed','primal_enriched','enriched_atoms','enriched_final')
RECEIPTS = ('range_two_replay','range_two_family_limit_replay','strict_family_separation',
            'two_spectator_closure','pair_transfer_overlap','one_spectator_overlap',
            'spin_overlap','coherent_overlap','charge_markov_extension','three_spectator_overlap')


def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p): return str(p.relative_to(ROOT))
def write(p,d): p.write_text(json.dumps(d,indent=2)+'\n')


def main():
    prior = read(BASE/'summary.json')
    # Existing reports and proofs remain intact; the previous full regression
    # applies to unchanged production sources, not to these added helper tests.
    prior_provenance = read(BASE/'provenance.json')['sha256']
    for name,h in prior_provenance.items():
        if sha(ROOT/name)!=h: raise ValueError('Changed preceding provenance: '+name)
    prior_entries=0
    for name in prior['current_receipts']:
        d=read(ROOT/name)
        if d.get('accepted') is not True: raise ValueError('Unaccepted preceding receipt')
        for source,h in d['source_sha256'].items():
            if sha(ROOT/source)!=h: raise ValueError('Stale preceding receipt source')
            prior_entries+=1
    if '7 passed' not in (BASE/'refinement_validation.log').read_text():
        raise ValueError('Seven focused refinement tests required')
    selected=BASE/'W_plus_1/enriched_final'
    for name in RECEIPTS:
        if read(selected/(name+'.json')).get('accepted') is not True:
            raise ValueError('Missing selected exact replay: '+name)
    paths=[BASE/'W_zero/final/three_spectator_overlap.json']
    files={Path(__file__).resolve(),BASE/'summary.json',BASE/'provenance.json'}
    for name in DIRECTORIES:
        directory=BASE/'W_plus_1'/name
        files.update(p for p in directory.rglob('*') if p.is_file())
        for receipt in RECEIPTS:
            p=directory/(receipt+'.json')
            if p.exists(): paths.append(p)
    entries=0
    for p in paths:
        d=read(p)
        if d.get('accepted') is not True: raise ValueError('Unaccepted new proof receipt')
        for source,h in d['source_sha256'].items():
            if sha(ROOT/source)!=h: raise ValueError('Stale new proof source: '+source)
            entries+=1;files.add(ROOT/source)
    energy=read(selected/'range_two_replay.json')
    family=read(selected/'range_two_family_limit_replay.json')
    old=prior['cases']['W_plus_1']
    lo=F(family['accepted_periodic_lower']);cap=F(family['periodic_family_upper']);gap=cap-lo
    if lo!=F(old['periodic_lower']) or gap!=F(family['family_gap']) or not 0<=gap<F(old['family_gap']):
        raise ValueError('Refined interval does not improve the preceding exact gap')
    if (selected/'profile_joint_r1_2_certificate.json').read_bytes()!=(BASE/'W_plus_1/final/profile_joint_r1_2_certificate.json').read_bytes():
        raise ValueError('Selected physical energy recipe changed')
    if any(map(F,read(selected/'two_spectator_closure.json')['two_spectator_moments'].values())):
        raise ValueError('Two-spectator closure failed')
    for name in ('pair_transfer_overlap','spin_overlap'):
        if any(F(x['exact_moment']) for x in read(selected/(name+'.json'))['separators']):
            raise ValueError('Preceding closure failed')
    if any(map(F,read(selected/'one_spectator_overlap.json')['spectator_moments'].values())) or F(read(selected/'coherent_overlap.json')['exact_moment']):
        raise ValueError('Spectator or coherent closure failed')
    probes={}
    for case,directory in [('W_zero',BASE/'W_zero/final'),('W_plus_1',selected)]:
        d=read(directory/'three_spectator_overlap.json')
        if d['directions']!=18 or d['violated_directions']!=18 or d['new_independent_directions']!=18:
            raise ValueError('Expected exact independent obstructions missing')
        largest=max(d['moments'],key=lambda x:abs(F(x['exact_moment'])))
        probes[case]={'receipt':rel(directory/'three_spectator_overlap.json'),
                      'directions':18,'violations':18,'old_rank':d['old_rank']['rank'],
                      'combined_rank':d['combined_rank']['rank'],'largest':largest}
    construction_entries=0
    for name in DIRECTORIES:
        for p in (BASE/'W_plus_1'/name).rglob('*.json'):
            d=read(p)
            if not isinstance(d,dict):continue
            for source,h in d.get('source_sha256',{}).items():
                if sha(ROOT/source)!=h:raise ValueError('Stale construction provenance: '+source)
                construction_entries+=1;files.add(ROOT/source)
    discovery=BASE.parent/'discovery'
    files.update(discovery.glob('two_spectator*.py'))
    files.update(discovery/n for n in ('three_spectator_overlap_probe.py','one_column_physical_completion.py'))
    files.update(ROOT/'tests'/n for n in ('test_marginal_three_spectator_probe.py','test_marginal_one_column_completion.py'))
    files.update(BASE/n for n in ('refinement_validation.log','three_spectator_validation.log','physical_completion_validation.log'))
    files.update([BASE/'W_plus_1/primal_subspace.log',BASE/'W_plus_1/primal_enriched.log',BASE/'W_zero/final/three_spectator_probe.log'])
    files.update(paths)
    report=ROOT/'research/marginal_two_spectator_refinement.md';files.add(report)
    diagnostics={}
    for name in DIRECTORIES:
        directory=BASE/'W_plus_1'/name
        for filename in ('diagonal_family_limit_diagnostic.json','subspace_diagnostic.json','pruning_diagnostic.json'):
            p=directory/filename
            if p.exists():
                d=read(p)
                diagnostics[rel(p)]={k:v for k,v in d.items() if k not in ('source_sha256','rational_family_upper','retained_indices','active_sectors')}
    summary={'accepted':True,'selected_directory':rel(selected),'periodic_lower':str(lo),
             'periodic_family_ceiling':str(cap),'family_gap':str(gap),'family_gap_float':float(gap),
             'previous_family_gap':old['family_gap'],'fraction_gap_closed':str(1-gap/F(old['family_gap'])),
             'fraction_gap_closed_float':float(1-gap/F(old['family_gap'])),
             'mixture_sources':family['family_replay']['nearest_constraint_replay']['mixture_sources'],
             'open_lower_per_site':energy['lower_per_site'],'physical_upper_per_site':energy['upper_per_site'],
             'three_spectator_obstructions':probes,'current_receipts':list(map(rel,paths)),
             'current_receipts_audited':len(paths),'current_source_hash_entries_checked':entries,
             'prior_current_receipts_reaudited':len(prior['current_receipts']),
             'prior_current_source_hash_entries_rechecked':prior_entries,
             'preceding_provenance_files_unchanged':len(prior_provenance),
             'construction_source_hash_entries_checked':construction_entries,
             'validation':{'unchanged_production_full_tests':945,'unchanged_production_full_subtests':102,
                           'additional_focused_tests':7,'production_changed':False,
                           'scope':'Previous full suite applies to unchanged production sources. Seven new discovery/completion tests passed separately; no new full-suite claim.'},
             'diagnostics':diagnostics,'report':rel(report),'provenance':rel(BASE/'refinement_provenance.json'),
             'remaining':['W1 finite gap remains nonzero; no numerical-limit convergence or exact attainment claim.',
                          'All18 new independent three-spectator conditions fail in both selected mixtures; stronger energy integration remains open.',
                          'No new coupling-transfer run in this refinement. Previous exact transfer preserved.',
                          'General quantum representability, generic molecular/long-range/higher-dimensional transfer and requested-accuracy scalability remain unproved. Goal active.']}
    write(BASE/'refinement_provenance.json',{'scope':'Current new proof and construction provenance. Numerical proposals remain nonaccepting. Earlier381-file provenance independently checked unchanged.',
                                          'sha256':{rel(p):sha(p) for p in sorted(files)}})
    write(BASE/'refinement_summary.json',summary)
    central=ROOT/'results/marginal_final_validation.json';d=read(central)
    d['two_spectator_refinement']=summary
    d['latest_validation_scope']='W1 family cap refined by candidate-feasible pricing and coherent subspace atoms, with exact physical replay and unchanged lower energy certificate. All18 independent three-spectator conditions violate both selected mixtures. Production sources unchanged from945 tests+102subtests;7 new focused tests pass. Numerical limit, stronger energy integration and general representability remain open.'
    write(central,d)
    print(json.dumps({'accepted':True,'family_gap':float(gap),'gap_reduction_fraction':summary['fraction_gap_closed_float'],
                      'new_receipts':len(paths),'new_source_entries':entries,'construction_entries':construction_entries,
                      'provenance_files':len(files),'prior_provenance_files_unchanged':len(prior_provenance)}))


if __name__=='__main__':main()
