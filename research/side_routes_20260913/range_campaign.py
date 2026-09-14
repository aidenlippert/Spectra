"""Bounded numerical search followed by exact rational training and transfer checks."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import platform
import socket
import time

import numpy as np
import scipy
from scipy.optimize import differential_evolution, minimize

from research.side_routes_20260913.chain_campaign import case
from research.side_routes_20260913.finite_range import bounds
from research.side_routes_20260913.range_discovery import NumericalChain, amplitude_from_weights, rational_amplitude, site_features


def objective_value(lower, upper, objective):
    if objective == 'width': return upper-lower
    if objective == 'energy': return upper
    if objective == 'lower': return -lower
    raise ValueError('Unknown search objective')


def search(task):
    name, radius, objective, seed = task
    objective_value(0,0,objective)
    started = time.monotonic()
    model = case(name, 12)
    kernel = NumericalChain(model, radius)
    compiled = time.monotonic()
    features = np.array(site_features(model))
    active = list(range(radius)) + [radius+j for j in range(3) if np.ptp(features[:,j]) > 0]
    initial = np.zeros(radius+3)
    if radius:
        initial[0] = np.log({'repulsive':.6,'site_disorder':2/3,'bond_disorder':.9}[name])
    best = [np.inf, initial.copy()]
    calls = 0

    def evaluate(short):
        nonlocal calls
        theta = np.zeros(radius+3)
        theta[active] = short
        lo, hi = kernel.interval(theta)
        score = objective_value(lo,hi,objective)
        calls += 1
        if score < best[0]: best[:] = score, theta.copy()
        return score

    evaluate(initial[active])
    de = differential_evolution(evaluate, [(-1.5,1.5)]*len(active), seed=seed,
                                popsize=6, maxiter=30, tol=1e-8, atol=0,
                                polish=False, x0=initial[active], workers=1)
    local = minimize(evaluate, best[1][active], method='Powell',
                     bounds=[(-1.5,1.5)]*len(active),
                     options={'maxfev':600,'xtol':1e-5,'ftol':1e-7})
    searched = time.monotonic()
    candidates = []
    for label, theta in (('initial',initial),('optimized',best[1])):
        amplitude, weights = rational_amplitude(model,radius,theta)
        claim = bounds(model,amplitude)
        proposed = kernel.interval(theta)
        rounded = kernel.interval(np.log([float(F(w)) for w in weights]))
        if max(abs(rounded[0]-claim['lower_float']),abs(rounded[1]-claim['upper_float'])) > 1e-8:
            raise AssertionError('Numerical/exact training disagreement after rationalization')
        candidates.append({'candidate':label,'weights':weights,'model':model,'amplitude':amplitude,'claim':claim,
                           'continuous_lower':proposed[0],'continuous_upper':proposed[1],
                           'rationalization_width_change':claim['width_float']-(proposed[1]-proposed[0])})
    winner = min(candidates,key=lambda p:objective_value(F(p['claim']['lower']),F(p['claim']['upper']),objective))
    return {'case':name,'radius':radius,'objective':objective,'seed':seed,'training_modes':12,
            'active_parameters':len(active),'numerical_objective_calls':calls,'compile_seconds':compiled-started,
            'numerical_search_seconds':searched-compiled,'wall_seconds':time.monotonic()-started,
            'de_reported_success':bool(de.success),'de_message':str(de.message),
            'local_reported_success':bool(local.success),'local_message':str(local.message),
            'winner':winner,'candidates':candidates,
            'scope':'Budgeted heuristic search. Optimizer completion is not proof of a global optimum.'}


def transfer(task):
    name, objective, radius, m, weights, origin = task
    model = case(name,m)
    amplitude = amplitude_from_weights(model,radius,weights)
    return {'case':name,'objective':objective,'radius':radius,'modes':m,'training_modes':12,
            'weights':weights,'origin':origin,'model':model,'amplitude':amplitude,'claim':bounds(model,amplitude)}


def run(out, workers, objectives=('width','energy')):
    out.mkdir(parents=True,exist_ok=False)
    (out/'searches').mkdir(); (out/'witnesses').mkdir()
    start = time.monotonic()
    names = ('repulsive','bond_disorder','site_disorder')
    tasks = [(name,r,obj,seed) for name in names for r in range(4)
             for obj in objectives for seed in (17,43)]
    searches = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending = {pool.submit(search,t):t for t in tasks}
        for future in as_completed(pending):
            record = future.result()
            searches.append(record)
            name,r,obj,seed = pending[future]
            (out/'searches'/f'{name}_{obj}_r{r}_seed{seed}.json').write_text(json.dumps(record,indent=2)+'\n')
            print('SEARCH',name,obj,r,seed,record['winner']['claim']['width_float'],record['numerical_objective_calls'],flush=True)
    search_finished = time.monotonic()
    selected = []
    transfer_tasks = []
    for name in names:
        for obj in objectives:
            for radius in range(4):
                eligible = [s for s in searches if s['case']==name and s['objective']==obj and s['radius']<=radius]
                winner = min(eligible,key=lambda s:objective_value(F(s['winner']['claim']['lower']),F(s['winner']['claim']['upper']),obj))
                r = winner['radius']; raw = winner['winner']['weights']
                weights = raw[:r]+['1']*(radius-r)+raw[r:]
                origin = {'search_radius':r,'seed':winner['seed'],'candidate':winner['winner']['candidate']}
                selected.append({'case':name,'objective':obj,'radius':radius,'weights':weights,'origin':origin})
                for m in (8,12,16,32,64):
                    transfer_tasks.append((name,obj,radius,m,weights,origin))
    witnesses = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending = {pool.submit(transfer,t):t for t in transfer_tasks}
        for future in as_completed(pending):
            result = future.result(); witnesses.append(result)
            name,obj,r,m,*_ = pending[future]
            (out/'witnesses'/f'{name}_{obj}_r{r}_m{m}.json').write_text(json.dumps(result,indent=2)+'\n')
            print('EXACT',name,obj,r,m,result['claim']['width_float'],result['claim']['wall_seconds'],flush=True)
    searches.sort(key=lambda s:(s['case'],s['objective'],s['radius'],s['seed']))
    witnesses.sort(key=lambda s:(s['case'],s['objective'],s['radius'],s['modes']))
    receipt = {'host':socket.gethostname(),'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,
               'workers':workers,'searches':len(searches),'exact_transfer_witnesses':len(witnesses),
               'numerical_objective_calls':sum(s['numerical_objective_calls'] for s in searches),
               'search_phase_wall_seconds':search_finished-start,'total_wall_seconds':time.monotonic()-start,
               'sum_search_job_wall_seconds':sum(s['wall_seconds'] for s in searches),
               'sum_transfer_exact_wall_seconds':sum(s['claim']['wall_seconds'] for s in witnesses),
               'training_modes':12,'parameter_quantization_denominator':65536,
               'selection':'Exact training objective across both seeds and all smaller radii, embedded with unit pair weights. Parameters then fixed for every transfer size.',
               'selected':selected,'ladder':[{k:s[k] for k in ('case','objective','radius','modes','origin')}|s['claim'] for s in witnesses],
               'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')},
               'scope':'Three abstract chain scenarios; heuristic discovery with finite search budgets. All accepted reported intervals use exact rational arithmetic.'}
    (out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('DONE',receipt['total_wall_seconds'],receipt['numerical_objective_calls'],flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--objectives',nargs='+',choices=('width','energy','lower'),default=['width','energy'])
    args = parser.parse_args()
    if not 1 <= args.workers <= 8: raise ValueError('Worker budget must be 1 through 8')
    run(args.out,args.workers,args.objectives)
