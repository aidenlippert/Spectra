"""Collect accepted completion results while retaining earlier checkpoints."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8/spectator_hopping'
CASES = ('W_zero', 'W_plus_1')
NAMES = ('range_two_replay','range_two_family_limit_replay','ceiling_improvement',
         'pair_transfer_overlap','one_spectator_overlap','spin_overlap',
         'coherent_overlap','charge_markov_extension')


def read(path):
    return json.loads(path.read_text())


def relative(path):
    return str(path.relative_to(ROOT))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data,indent=2)+'\n')


def main():
    earlier = read(BASE/'combined_summary.json')
    prior = read(BASE/'limit_refinement_summary.json')
    old_paths = earlier['receipts'] + prior['new_receipts'] + prior['supplemental_receipts']
    new_paths = [relative(BASE/c/'completed_basis'/(name+'.json')) for c in CASES for name in NAMES]
    count = 0
    for name in old_paths + new_paths:
        data = read(ROOT/name)
        assert data['accepted'], name
        for source, expected in data['source_sha256'].items():
            assert sha(ROOT/source) == expected, (name, source)
            count += 1
    cases = {}
    for case in CASES:
        folder = BASE/case/'completed_basis'
        energy = read(folder/'range_two_replay.json')
        family = read(folder/'range_two_family_limit_replay.json')
        diagnostic = read(folder/'basis_completion_diagnostic.json')
        improvement = read(folder/'ceiling_improvement.json')
        proposal = read(folder/'diagonal_family_limit_proposal.json')
        assert diagnostic['proposal_written'] and not diagnostic['accepted']
        assert len(proposal['mixture']) == 85
        for source, expected in diagnostic['source_sha256'].items():
            assert sha(ROOT/source) == expected, source
        assert energy['lower_per_site'] == prior['cases'][case]['lower_per_site']
        assert energy['upper_per_site'] == prior['cases'][case]['upper_per_site']
        pair = read(folder/'pair_transfer_overlap.json')
        assert all(item['violated'] for item in pair['separators'])
        assert all(F(v)==0 for v in read(folder/'one_spectator_overlap.json')['spectator_moments'].values())
        cases[case] = {
            'directory':relative(folder), 'target':energy['target'],
            'lower_per_site':energy['lower_per_site'], 'upper_per_site':energy['upper_per_site'],
            'periodic_lower':family['accepted_periodic_lower'],
            'periodic_family_ceiling':family['periodic_family_upper'],
            'family_ceiling_float':float(F(family['periodic_family_upper'])),
            'family_gap':family['family_gap'], 'family_gap_float':family['family_gap_float'],
            'exact_ceiling_improvement':improvement['exact_improvement'],
            'ceiling_improvement_float':improvement['improvement_float'],
            'gap_reduction_fraction':improvement['gap_reduction_fraction'],
            'mixture_sources':85,'selected_columns':diagnostic['selected_columns'],
            'completions_tried':diagnostic['completions_tried'],
            'nonnegative_completions':diagnostic['nonnegative_completions'],
            'added_old_columns':diagnostic['added_old_columns'],
            'maximum_weight_string_length':max(len(item['weight']) for item in proposal['mixture']),
            'pair_transfer_separators':pair['separators'],
        }
    test_paths = [ROOT/'results/marginal_graded_hubbard8/discovery/fraction_free_completion.py',
                  ROOT/'tests/test_marginal_fraction_free_completion.py',
                  BASE/'basis_completion_focused_validation.log']
    assert '9 passed in 0.03s' in test_paths[-1].read_text()
    focused = {
        'command':'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q tests/test_marginal_fraction_free_completion.py',
        'tests_passed':9,'elapsed_seconds':.03,'exit_code':0,
        'source_sha256':{relative(p):sha(p) for p in test_paths},
        'scope':'New integer elimination helper: original-equation reconstruction, row swaps, negative determinant, inconsistent rectangular systems, explicit residual completion and refusal cases. No full-suite rerun.',
    }
    write(BASE/'basis_completion_focused_validation.json',focused)
    files = set(test_paths)
    for case in CASES:
        folder = BASE/case/'completed_basis'
        files.update(p for p in folder.iterdir() if p.is_file())
        files.add(BASE/case/'scaled_basis/candidate_ledger.json.gz')
    for name in ('fraction_free_completion.py','spectator_basis_completion.py',
                 'spectator_basis_completion_anchors.py','spectator_ceiling_improvement.py'):
        files.add(ROOT/'results/marginal_graded_hubbard8/discovery'/name)
    write(BASE/'basis_completion_provenance.json',{
        'scope':'Artifact provenance, not an independent physical acceptance gate. W0 uses spectator_basis_completion.py; W1 uses the determinant-anchor variant. Both use saved ledgers and the prior polished85-source pool. No spectral sampling or pricing was rerun.',
        'sha256':{relative(p):sha(p) for p in sorted(files)},
    })
    summary = {
        'accepted':True,'cases':cases,'new_receipts':new_paths,
        'all_receipts_audited':old_paths+new_paths,
        'receipts_audited':len(old_paths)+len(new_paths),
        'source_hash_entries_checked':count,'focused_validation':focused,
        'production_verifier_changed':False,
        'report':'research/marginal_spectator_basis_completion.md',
        'provenance':relative(BASE/'basis_completion_provenance.json'),
        'scope':'Exact85-source ceilings close to the previously inconsistent numerical proposals, with strict improvement at unchanged accepted lower certificates. Fixed current family only; no physical energy upper or optimum claim.',
        'remaining':[
            'The unrestricted family optimum remains unresolved; current gaps are about5.73e-5 and4.64e-4.',
            'All four independent pair-transfer violations persist in both new mixtures; energy integration remains undone.',
            'W=+/-0.1 retain their older polished checkpoints. No new transfer beyond these matched targets was established.',
            'General quantum representability and requested-accuracy scalability remain unproved. Goal active.',
        ],
    }
    write(BASE/'basis_completion_summary.json',summary)
    status_path = ROOT/'results/marginal_final_validation.json'
    status = read(status_path)
    status['spectator_basis_completion'] = summary
    status['latest_validation_scope'] = 'Nine focused integer-completion tests and16fresh accepting receipts;81historical/current receipts hash-audited. Production verifiers unchanged. Prior847-test/102-subtest full regression was not rerun.'
    write(status_path,status)
    print(json.dumps({'receipts':len(old_paths)+len(new_paths),'hash_entries':count,'focused_tests':9,
                      'cases':{c:{'ceiling':d['family_ceiling_float'],'gap':d['family_gap_float']} for c,d in cases.items()}}))


if __name__ == '__main__':
    main()
