"""Replay combined bounds and compare them with independent small references."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import time

from research.side_routes_20260913.combine_bounds import combine


def analyze_group(task):
    name,radius,m,candidates,reference = task
    paths,payloads = zip(*candidates)
    result = combine(payloads)
    low,high = result['lower_source_index'],result['upper_source_index']
    chosen = sorted({low,high})
    result.update({'case':name,'maximum_radius':radius,'modes':m,
                   'lower_source':paths[low],'upper_source':paths[high],
                   'lower_objective':payloads[low]['objective'],'upper_objective':payloads[high]['objective'],
                   'lower_radius':payloads[low]['radius'],'upper_radius':payloads[high]['radius'],
                   'evaluated_candidate_count':len(payloads),
                   'selected_amplitude_json_bytes':sum(len(json.dumps(payloads[i]['amplitude']).encode()) for i in chosen),
                   'selected_payload_sha256':{paths[i]:hashlib.sha256(json.dumps(payloads[i],sort_keys=True).encode()).hexdigest() for i in chosen},
                   'best_single_witness_width':str(min(F(p['claim']['width']) for p in payloads))})
    if reference is not None:
        if reference['model'] != payloads[0]['model']:
            raise ValueError('Reference Hamiltonian differs')
        energy = reference['energy_float']
        lower_slack,upper_slack = energy-result['lower_float'],result['upper_float']-energy
        if min(lower_slack,upper_slack) < -1e-8:
            raise AssertionError('Numerical reference lies outside the exact interval')
        result['numerical_reference'] = {'energy':energy,'residual_norm':reference['residual_norm'],
                                          'lower_slack':lower_slack,'upper_slack':upper_slack,
                                          'lower_fraction_of_width':lower_slack/result['width_float']}
    return result


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--campaign',action='append',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workers',type=int,default=4)
    args=p.parse_args();started=time.monotonic()
    if not 1<=args.workers<=8: raise ValueError('Invalid worker budget')
    args.out.mkdir(parents=True,exist_ok=False)
    candidates=[]
    for directory in args.campaign:
        candidates.extend((str(path),json.loads(path.read_text())) for path in sorted((directory/'witnesses').glob('*.json')))
    references=json.loads(args.reference.read_text())['references']
    tasks=[]
    for name in ('repulsive','bond_disorder','site_disorder'):
        for radius in range(4):
            for m in (8,12,16,32,64):
                choices=[item for item in candidates if item[1]['case']==name and item[1]['radius']<=radius and item[1]['modes']==m]
                reference=next((r for r in references if r['case']==name and r['modes']==m),None)
                tasks.append((name,radius,m,choices,reference))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows=list(pool.map(analyze_group,tasks))
    for row in rows:
        (args.out/f"{row['case']}_r{row['maximum_radius']}_m{row['modes']}.json").write_text(json.dumps(row,indent=2)+'\n')
    receipt={'combined_certificates':len(rows),'rows':rows,'wall_seconds':time.monotonic()-started,
             'all_selected_certificates_replayed':True,'reference_sha256':hashlib.sha256(args.reference.read_bytes()).hexdigest(),
             'scope':'Max lower and min upper across evaluated objectives and radii, for identical models. Parameters were fitted at M12; selection of bounds uses their exact values at each size.'}
    (args.out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'combined_certificates':len(rows),'wall_seconds':receipt['wall_seconds']}))
