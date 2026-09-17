"""Frozen strong-control transfer experiment, using only compact supplied inputs.

The selected population contrast is between the most and least occupied
spatial orbitals. This changes the observable between molecules and is a
declared construction rule, not a matched physical-control comparison.
"""
import argparse,json,time,sys
from fractions import Fraction as F
from pathlib import Path
from research.correlated_pair_20260913.mps_exact import State
from .collective_control import construct,replay


def choose_pair(data,cert):
    state=State(data,cert)
    occupations=[]
    for p in range(0,data['modes'],2):
        occupations.append(state.expectation({((1,p),(0,p)):F(1),((1,p+1),(0,p+1)):F(1)}))
    high=max(range(len(occupations)),key=lambda j:(occupations[j],-j))
    low=min((j for j in range(len(occupations)) if j!=high),key=lambda j:(occupations[j],j))
    return (2*high,2*low),occupations,state.stats


def run(output):
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    start=time.monotonic();rows=[]
    for case in ('h6_asymmetric','water_asymmetric'):
        path=Path('results/transfer_solver_20260915/cases')/case
        data=json.loads((path/'fixture.json').read_text());cert=json.loads((path/'mps/state.json').read_text())
        st=time.monotonic();pair,occ,stats=choose_pair(data,cert)
        artifact=construct(data,cert,max_amplitude=F(1024),pair=pair)
        construction=time.monotonic()-st
        (output/(case+'.json')).write_text(json.dumps(artifact,indent=2)+'\n')
        st=time.monotonic();pair2,occ2,_=choose_pair(data,cert)
        if pair2!=pair or occ2!=occ:raise AssertionError('Selection did not reproduce')
        fresh=replay(data,cert,artifact,F(1024),expected_pair=pair2)
        row={'case':case,'fixture':str(path/'fixture.json'),'state':str(path/'mps/state.json'),
             'control_pair_even_spin_indices':list(pair),'spatial_occupations':[str(x) for x in occ],
             'construction_and_selection_seconds':construction,'fresh_replay_and_selection_seconds':time.monotonic()-st,
             'selection_stats':stats,'certificate_stats':artifact['stats'],
             'selected_protocol':fresh['attempts'][-1],'accepted':fresh['accepted']}
        rows.append(row)
        print(json.dumps({'case':case,'accepted':row['accepted'],'pair':pair,'construction_seconds':construction,
                          'amplitude':fresh['attempts'][-1]['v_Ha'],'D_interval':fresh['attempts'][-1]['D_interval_float']}),flush=True)
    if any(x in sys.modules for x in ('numpy','scipy','quimb','pyscf')):raise AssertionError('Exact path imported numeric library')
    result={'kind':'collective_strong_control_transfer_v1','rule':'Max/min exact spatial occupation with smallest-index tie break; theta=25/8; v doubles from 1/2 to at most 1024 Ha; target D<=-3/5; original short-case uncertainty budgets.',
            'cases':rows,'seconds':time.monotonic()-start,'all_accepted':all(r['accepted'] for r in rows),
            'scope':'Conditional strong mathematical control only. Observables selected by the frozen occupancy rule. Supplied model and initial-state discovery inherited. Not the original low-amplitude problem or laboratory transfer.'}
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.output)
