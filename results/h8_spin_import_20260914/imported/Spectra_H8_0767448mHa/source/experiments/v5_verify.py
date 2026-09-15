"""Verify saved V5 traces, mathematical arithmetic and resource summaries.

Trace replay is a reproducibility check, not an independent proof of the source
law. Evidence reconstruction and counter arithmetic also operate on saved data.
"""
from pathlib import Path
from fractions import Fraction as F
from dataclasses import asdict
from hashlib import sha256
import json
import re
from experiments.v5_core import replay_records, PREPARATION_TIME
from experiments.v5_run import METHODS, run_world, world, serialize
from experiments.v5_complementarity import certificate, run as run_retention

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {'cumulative':21,'frozen_after_first':37,'scratch':56,'exact_lookup':56,
            'structured_retrieval':21,'stateful_bayes':21,'fixed_one_read':20}


def without_timing(row):
    return {k:v for k,v in row.items() if k != 'cpu_seconds'}


def verify():
    saved=json.loads((ROOT/'results/v5/transfer_results.json').read_text())
    source_reads=tasks=artifacts=0
    for row in saved['world_results']:
        data=world(row['evaluator_seed'])
        for method in METHODS:
            result=row['methods'][method]
            replay=serialize(run_world(data,method))
            assert without_timing(result)==without_timing(replay), 'trace replay mismatch'
            assert result['total_scalar_reads']==EXPECTED[method]
            evidence=set()
            for source,ledger in zip(result['source_records'],result['source_ledgers']):
                n=len(source['records'])
                assert ledger['scalar_readouts']==ledger['binary_records']==ledger['reset_operations']==n
                assert ledger['source_preparations']==n and ledger['source_coordinate_preparations']==4*n
                assert ledger['relaxation_time_units']==n*PREPARATION_TIME
                assert not ledger['stationary_reference'] and not ledger['intentional_contract_violation']
                for rec in source['records']:
                    evidence.add((tuple(source['lineage'][:rec['prefix']]),rec['alpha'],
                                  tuple(rec['probe']),rec['outcome']))
                source_reads+=n
            for task in result['tasks']:
                assert task['correct']==(task['estimate']==task['evaluator_hidden_prefix'][-1])
                tasks+=1
                if 'learned' in task:
                    learned=task['learned']
                    rebuilt=replay_records(learned['records'])
                    assert serialize(asdict(rebuilt))==learned
                    assert F(learned['failure_bound'])<F(1,20)
                    for rec in learned['records']:
                        assert (tuple(rec['lineage']),rec['alpha'],tuple(rec['probe']),rec['outcome']) in evidence
                    artifacts+=1
            if 'initial_artifact' in result:
                assert serialize(asdict(replay_records(result['initial_artifact']['records'])))==result['initial_artifact']
    # Recompute summaries rather than trusting runner annotations.
    for method,summary in saved['summary'].items():
        rows=[r['methods'][method] for r in saved['world_results']]
        assert summary['source_reads']==sum(r['total_scalar_reads'] for r in rows)
        assert summary['relaxation_time_units']==summary['source_reads']*PREPARATION_TIME
        for stage in (2,3):
            target_rows=[t for r in rows for t in r['tasks'] if t['stage']==stage]
            assert summary[f'stage_{stage}_tasks']==len(target_rows)
            assert summary[f'stage_{stage}_errors']==sum(not t['correct'] for t in target_rows)
    bounds=saved['analytical_bounds']
    r0=F(123,496)-F(1,10**9)
    lower=F(9,10)*(1+2*r0)
    assert F(bounds['adaptive_expected_cold_reads_lower_bound'])==lower>F(5,4)
    assert F(bounds['net_expected_read_saving_per_four_target_batch_lower_bound'])==4*lower-5>0
    retention=json.loads((ROOT/'results/v5/complementarity_results.json').read_text())
    assert retention==json.loads(json.dumps(run_retention(retention['worlds'],retention['targets_per_world'])))
    assert retention['certificate']==json.loads(json.dumps(certificate()))
    assert F(retention['certificate']['acquired_complementarity_lower'])>F(77,1000)
    profile=json.loads((ROOT/'results/v5/computational_accounting.json').read_text())
    read_count=next(x['total_calls'] for x in profile['selected_function_counts']
                    if x['file']=='v5_core.py' and x['function']=='read')
    assert read_count==source_reads+retention['acquisition_reads']+retention['counterfactual_target_reads']
    assert all(m['max_record_coordinate_bits']<=8 for row in saved['world_results'] for m in row['methods'].values())
    log=(ROOT/'results/v5/test_log.txt').read_text()
    assert log.rstrip().endswith('OK')
    tests=int(re.search(r'Ran (\d+) tests',log).group(1))
    code=sorted(list(ROOT.glob('experiments/v5_*.py'))+list(ROOT.glob('tests/test_v5*.py')))
    result_paths=sorted((ROOT/'results/v5').glob('*results.json'))
    result_paths.append(ROOT/'results/v5/computational_accounting.json')
    documents=[ROOT/'research/v5/PROOF.md', ROOT/'research/v5/goal_audit.md', ROOT/'results/v5/REPORT.md']
    receipt={'test_count':tests,'tests_passed':True,'saved_transfer_source_reads_replayed':source_reads,
             'saved_target_decisions_checked':tasks,'learned_artifacts_rebuilt':artifacts,
             'transfer_worlds':len(saved['world_results']),'retention_worlds':retention['worlds'],
             'exact_positive_complementarity_lower':retention['certificate']['acquired_complementarity_lower'],
             'two_acquisition_inclusive_gains_positive':True,'structured_conventional_ties':True,
             'profile_source_calls_reconciled':read_count,
             'physical_experiments_performed':False,'unrestricted_scientific_capability_proved':False,
             'limits':['Source and reset availability are physical model assumptions.',
                       'The model grammar and experiment-compilation algebra are programmed.',
                       'Numerical trace replay does not prove the continuous stochastic law.',
                       'Complete hardware energy and general bit-complexity bounds are absent.'],
             'source_hashes':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in code},
             'document_hashes':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in documents},
             'result_hashes':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in result_paths}}
    (ROOT/'results/v5/verification_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    print(json.dumps({k:v for k,v in verify().items() if not k.endswith('hashes')},indent=2))
