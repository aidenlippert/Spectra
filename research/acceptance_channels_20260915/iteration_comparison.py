"""Read matched logs without treating a time-budget result as a family advantage."""
import hashlib
import json
from research.acceptance_channels_20260915.campaign import OUT,dump


def records(path):
    result={}
    for line in path.read_text().splitlines():
        try:row=json.loads(line)
        except ValueError:continue
        if type(row.get('iteration')) is int and 'unverified_width_mHa' in row:
            if row['iteration'] in result:raise ValueError('Repeated iteration in log')
            result[row['iteration']]=row
    if 0 not in result:raise ValueError('Initial checkpoint score missing')
    return result


def main():
    paths={name:OUT/'runs'/('h12_'+tag+'.log') for name,tag in
        (('control','matched_continuation'),('four','dense_t2_solve'))}
    logs={name:records(path) for name,path in paths.items()}
    if logs['control'][0]['unverified_width_mHa']!=logs['four'][0]['unverified_width_mHa']:
        raise ValueError('Initial checkpoint scores differ')
    rows=[]
    for i in sorted(set(logs['control'])&set(logs['four'])):
        control,four=logs['control'][i],logs['four'][i]
        rows.append({'iteration':i,'control_predicted_width_mHa':control['unverified_width_mHa'],
            'four_predicted_width_mHa':four['unverified_width_mHa'],
            'control_minus_four_mHa':control['unverified_width_mHa']-four['unverified_width_mHa'],
            'control_elapsed_seconds':control['seconds'],'four_elapsed_seconds':four['seconds']})
    discovered={name:json.loads((OUT/'h12'/tag/'discovery.json').read_text()) for name,tag in
        (('control','matched_continuation'),('four','dense_t2'))}
    result={'scope':'Read-only comparison of numerical iterates; not new exact certificates or a family obstruction.',
        'source_log_sha256':{name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in paths.items()},
        'common_logged_iterates':rows,'last_common_logged_iterate':rows[-1],
        'actual_completed_iterations':{name:d['iterations_completed'] for name,d in discovered.items()},
        'same_numerical_time_budget_seconds':300,
        'interpretation':'The four-direction run was worse at every positive common logged iteration. Its smaller final time-budget export does not isolate a benefit from added correlations; it completed more iterations. No reproducible throughput advantage has been established.'}
    dump(OUT/'iteration_comparison.json',result)
    print(json.dumps({key:result[key] for key in ('last_common_logged_iterate','actual_completed_iterations','interpretation')}),flush=True)


if __name__=='__main__':main()
