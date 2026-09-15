"""Collect exact interval receipts and complete cost ledgers for this pass."""
from fractions import Fraction as F
from pathlib import Path
import contextlib
import hashlib
import io
import json
import sys

from experiments.marginal_symbolic import decode
from research.certificate_scaling.cubic_interval_replay import run as replay_interval

ROOT=Path(__file__).resolve().parents[2]
RESULT=ROOT/'results/molecular_identity_20260913'


def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    experiments={}
    for name in ('h4_adaptive','h6_adaptive','h4_coupled','h6_coupled'):
        experiments[name]=read(RESULT/name/'experiment.json')
    if experiments['h4_coupled']['policy']!=experiments['h6_coupled']['policy']:
        raise ValueError('Coupled policy changed between H4 and H6')
    if experiments['h4_adaptive']['policy']!=experiments['h6_adaptive']['policy']:
        raise ValueError('Separate-block policy changed between H4 and H6')
    certificates=sorted(RESULT.glob('*/certificate.json'))+sorted(RESULT.glob('*/*/certificate.json'))
    rows=[]
    for certificate in certificates:
        cert=read(certificate);m=cert['modes'];label={8:'h4',12:'h6'}[m]
        fixture=ROOT/f'results/certificate_scaling/active_space_ladder/{label}/fixture.json'
        f=read(fixture)
        if cert['particles']!=f['particles'] or decode(cert['hamiltonian'],m,4)!=decode(f['hamiltonian'],m,4):
            raise ValueError('Certificate differs from frozen fixture')
        reference=ROOT/f'results/certificate_scaling/active_space_ladder_references_aligned/{label}/upper.json'
        interval_path=certificate.parent/'interval.json'
        if interval_path.exists():
            interval=read(interval_path)
            if interval['certificate_sha256']!=sha(certificate) or interval['reference_sha256']!=sha(reference):
                raise ValueError('Previously replayed input changed')
        else:
            with contextlib.redirect_stdout(io.StringIO()):
                interval=replay_interval(certificate,reference,interval_path)
        receipt=read(certificate.parent/'receipt.json')
        if F(receipt['exact']['lower'])!=F(interval['lower']):raise ValueError('Export/replay lower mismatch')
        row={'run':str(certificate.parent.relative_to(RESULT)),
             'lower':interval['lower'],'upper':interval['upper'],'width':interval['width'],
             'width_millihartree':1000*interval['width_float'],'passes_target':interval['passes_0_0016_Ha'],
             **{k:receipt[k] for k in ('polynomial_Gram_entries','gram_map_nonzeros','coefficient_rows',
                 'symbolic_word_pair_products','construction_seconds','solve_seconds','export_seconds','wall_seconds',
                 'status','certificate_bytes','factor_denominator_bits')},
             'factor_rows':receipt['exact']['factor_rows'],'factor_nonzeros':receipt['exact']['factor_nonzeros'],
             'max_factor_numerator_bits':max((abs(v).bit_length() for block in cert['blocks'] for factor in block['factor'] for v in factor),default=0),
             'residual_l1':receipt['exact']['residual_l1'],'residual_l1_float':receipt['exact']['residual_l1_float'],
             'interval_replay_seconds':interval['wall_seconds'],
             'upper_witness_states':interval['upper_replay']['witness_states'],
             'upper_word_state_checks':interval['upper_replay']['word_state_checks']}
        rows.append(row)
        print(json.dumps({k:row[k] for k in ('run','width_millihartree','passes_target','polynomial_Gram_entries')}),flush=True)
    costs={}
    for name,experiment in experiments.items():
        solved=[row for row in rows if row['run'].startswith(name+'/')]
        pricing=sum(item['pricing']['seconds'] for item in experiment['history'])
        counts=sum(item['pricing']['monomial_word_pair_products'] for item in experiment['history'])
        current={'solve_runs':len(solved),'pricing_seconds':pricing,'pricing_word_pair_products':counts,
                 'construction_seconds':sum(x['construction_seconds'] for x in solved),
                 'solve_seconds':sum(x['solve_seconds'] for x in solved),
                 'export_seconds':sum(x['export_seconds'] for x in solved),
                 'interval_replay_seconds':sum(x['interval_replay_seconds'] for x in solved),
                 'best_run':min(solved,key=lambda x:F(x['width']))['run']}
        if 'total_wall_seconds' in experiment:
            current['discovery_total_wall_seconds']=experiment['total_wall_seconds']
            if experiment.get('resumed_after_channel_fix'):
                current['timing_caveat']='Includes a pre-failure wall estimate recovered from timestamps; implementation-debugging idle time is excluded.'
                current['pre_failure_estimated_seconds']=experiment['prior_wall_seconds_estimated']
        else:
            current['discovery_total_wall_seconds']=experiment['additional_wall_seconds']+(experiment['seed_discovery_seconds'] if experiment['seed_reused'] else 0)
            current['reused_seed_seconds']=experiment['seed_discovery_seconds'] if experiment['seed_reused'] else 0
        current['export_scope']='Includes exact lower verification inside each solve.'
        costs[name]=current
    references={label:read(ROOT/f'results/certificate_scaling/active_space_ladder_references_aligned/{label}/receipt.json') for label in ('h4','h6')}
    historical={}
    for label in ('h4','h6'):
        p=ROOT/f'results/certificate_scaling/commutator_dictionary/{label}_creator_channels'
        historical[label]={'interval':read(p/'interval.json'),'receipt':read(p/'receipt.json')}
    result={'rows':rows,'policy_costs':costs,'reference_discovery':references,
            'historical_creator_baselines':historical,'both_transfer_policies_unchanged':True,
            'exact_intervals_accepted':len(rows),'scope':'Frozen rational Hamiltonians only. Timings are single campaign observations, with some concurrent runs; no controlled speedup claim.'}
    (RESULT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':run()
