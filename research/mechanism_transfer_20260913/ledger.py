"""Costs, canonical operator dimensions, precision and preservation receipt."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import time

from research.molecular_collective_20260913.core import extract,factor_operators
from research.mechanism_transfer_20260913.core import generator,combine

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'results/mechanism_transfer_20260913'


def modular_rank(rows):
    modulus=2**61-1
    if not rows:return 0
    a=[[int(F(x).numerator%modulus)*pow(F(x).denominator%modulus,-1,modulus)%modulus for x in row] for row in rows]
    pivot=0
    for j in range(len(a[0])):
        k=next((i for i in range(pivot,len(a)) if a[i][j]),None)
        if k is None:continue
        a[pivot],a[k]=a[k],a[pivot];inverse=pow(a[pivot][j],-1,modulus)
        a[pivot]=[(x*inverse)%modulus for x in a[pivot]]
        for i in range(pivot+1,len(a)):
            if a[i][j]:
                c=a[i][j];a[i]=[(x-c*y)%modulus for x,y in zip(a[i],a[pivot])]
        pivot+=1
        if pivot==len(a):break
    return pivot


def run():
    start=time.monotonic();prior=ROOT/'results/trace_pricing_20260913';mp=prior/'manifest.json'
    old=json.loads(mp.read_text())['files']
    changed=[n for n,v in old.items() if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=v['sha256']]
    if changed:raise ValueError('Frozen branch changed: '+repr(changed))
    watchdogs=json.loads((OUT/'watchdogs.json').read_text())+json.loads((OUT/'recovery_watchdogs.json').read_text())
    costs={'fixture_generation':0.,'discovery':0.}
    for row in watchdogs:costs[row['phase']]+=row['wall_seconds']
    cases=[]
    for directory in sorted((OUT/'campaign').iterdir()):
        if not directory.is_dir() or not (directory/'summary.json').exists():continue
        summary=json.loads((directory/'summary.json').read_text());rows=[]
        for trial in summary['trials']:
            row=dict(trial)
            if trial['accepted']:
                receipt=json.loads((directory/trial['arm']/'receipt.json').read_text())
                cert=json.loads((directory/trial['arm']/'certificate.json').read_text())
                row.update({k:receipt[k] for k in ('Gram_entries','anti_dimensions','certificate_bytes','constructed_map_nonzeros',
                    'construction_seconds','solve_seconds','export_seconds','accept_seconds','status')})
                row['expanded_certificate_bytes']=receipt['accepted']['expanded_certificate_bytes']
                blocks=cert['base_blocks']+cert['anti_blocks']
                factors=[x for b in blocks for v in b['factor'] for x in v if x]
                directions=[x for b in cert['anti_blocks'] for v in b['directions'] for x in v if x]
                row.update(factor_rows=sum(len(b['factor']) for b in blocks),factor_nonzeros=len(factors),
                    direction_nonzeros=len(directions),maximum_factor_integer_bits=max(abs(x).bit_length() for x in factors),
                    maximum_direction_integer_bits=max((abs(x).bit_length() for x in directions),default=0),
                    factor_denominator=cert['denominator'],direction_denominator=10**10,
                    exact_residual_penalty_mHa=float(1000*F(receipt['accepted']['retained']['residual_l1'])))
                factor_ranks=[modular_rank(b['factor']) for b in blocks]
                row['factor_rank_sum']=sum(factor_ranks)
                row['factor_rank_is_exact']=all(rank==len(b['factor']) for rank,b in zip(factor_ranks,blocks))
            rows.append(row)
        data=json.loads((directory/'fixture.json').read_text());tail=json.loads((directory/'tail.json').read_text())
        cert_path=directory/'coupled/certificate.json';ranks=[]
        if cert_path.exists():
            cert=json.loads(cert_path.read_text());p=extract(data,tail['center_number']);patterns=[q for _,q in factor_operators(p,tail)]
            for block in cert['anti_blocks']:
                base=[generator(patterns,ref,p['modes']) for ref in block['generators']]
                polys=[combine(base,v,block['direction_denominator']) for v in block['directions']]
                words=sorted(set(w for poly in polys for w in poly));matrix=[[poly.get(w,F(0)) for w in words] for poly in polys]
                rank=modular_rank(matrix)
                if rank!=len(polys):raise ValueError('Full canonical operator independence was not established')
                ranks.append({'dimension':len(polys),'modular_rank':rank,'canonical_monomials':len(words)})
        total=sum(x['wall_seconds'] for x in watchdogs if x['phase']=='discovery' and x['case']==summary['case'])
        cases.append({'case':summary['case'],'discovery_seconds_including_failed_attempts':total,
            'model_construction':summary['model_construction'],'rule_construction':summary['rule_construction'],
            'tail':summary['tail'],'canonical_operator_rank_checks':ranks,'trials':rows})
    analysis=json.loads((OUT/'training_analysis.json').read_text());old_cost=json.loads((prior/'cost_ledger.json').read_text())
    result={'preserved_prior_files':len(old),'preserved_manifest_sha256':hashlib.sha256(mp.read_bytes()).hexdigest(),
        'watchdogs':watchdogs,'new_case_costs_seconds':costs,'training_analysis_seconds':analysis['wall_seconds'],
        'inherited_training_certificate_discovery_seconds':old_cost['new_trace_cumulative_causal_discovery_seconds'],
        'new_discovery_plus_training_analysis_seconds':costs['discovery']+analysis['wall_seconds'],
        'cumulative_causal_discovery_seconds':old_cost['new_trace_cumulative_causal_discovery_seconds']+costs['discovery']+analysis['wall_seconds'],
        'cases':cases,'modular_rank_method':'Rational coefficients reduced modulo 2^61-1 using invertible denominators and pivots. Full row rank modulo this integer proves full row rank over Q; no primality assumption is needed for this conclusion.',
        'rank_and_ledger_seconds':time.monotonic()-start,
        'scope':'Discovery includes all three controls, failed setup attempts and corrected cases. Prior full-family diagnostics remain in the frozen ledger. Timed cases ran sequentially; verification follows discovery.'}
    if (OUT/'fresh_replay.json').exists():result['fresh_interval_replay_seconds']=json.loads((OUT/'fresh_replay.json').read_text())['wall_seconds']
    result['identity_check_seconds']=json.loads((OUT/'identity_check.json').read_text())['wall_seconds']
    result['verification_overlap']='The 0.657-second exact identity check overlapped fresh interval replay after all discovery. Do not add those wall times as calendar elapsed.'
    (OUT/'cost_ledger.json').write_text(json.dumps(result,indent=2)+'\n')
    files={ROOT/n for n in old};files.add(mp)
    files.update(p for p in (ROOT/'research/mechanism_transfer_20260913').rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    files.update(p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json')
    manifest={str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)}
    (OUT/'manifest.json').write_text(json.dumps({'files':manifest},indent=2)+'\n')


if __name__=='__main__':run()
