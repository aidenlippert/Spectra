"""Complete a deficient numerical basis using accepted physical source columns.

No spectral search. Candidate selection is numerical; elimination, completion,
positivity and every retained equation are exact. Production replay is separate.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import lcm
import argparse
import gzip
import hashlib
import json
import sys
import time
import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from signed_charge_numeric import prepare, BASE
from fraction_free_completion import reduce_basis
from spectator_family_numeric import solve_basis
from experiments.marginal_hopping_telescope import actions, projected_matrix
from experiments.marginal_spin_telescope import LABELS, actions as spin_actions
from experiments.marginal_spectator_hopping import LABELS as SPECTATORS, actions as spectator_actions
from experiments.marginal_joint_family_limit import _physical_vector
from experiments.marginal_signed_charge_telescope import PATTERNS,local_value as signed_value
from experiments.marginal_quadratic_charge_telescope import local_value as quadratic_value
from experiments.marginal_charge_square_pairs import local_value as square_value
from experiments.marginal_charge_indicator_telescope import local_value as indicator_value
from experiments.marginal_range_two_density import diagonal_value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--ledger', type=Path, required=True)
    parser.add_argument('--accepted-family', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    out = args.directory.resolve()
    cp = out / 'profile_joint_r1_2_certificate.json'
    c = json.loads(cp.read_text())
    old_path = args.accepted_family.resolve()
    old = json.loads(old_path.read_text())
    old_receipt_path = old_path.with_name('range_two_family_limit_replay.json')
    old_receipt = json.loads(old_receipt_path.read_text())
    if not old_receipt['accepted'] or old_receipt['source_sha256'][str(old_path.relative_to(ROOT))] != hashlib.sha256(old_path.read_bytes()).hexdigest():
        raise ValueError('Accepted old physical source set required')
    if len(old['mixture']) != 85 or any(F(item['weight']) <= 0 for item in old['mixture']):
        raise ValueError('Strictly positive85-source completion pool required')
    candidate_path = BASE / 'joint_projector/signed_density/symmetry_diagonal_candidates.json'
    candidates = json.loads(candidate_path.read_text())['candidates']
    shapes = [{int(s): v for s, v in candidates[i]['diagonal'].items()} for i in [1,2,4,5,6,7,31,11,20]]
    x, _, mats, half, hn, charged, cn, ratio = prepare(c, shapes)
    th = F(c['projector_sum_ceiling']) / c['windows']
    tj = F(c['joint']['projector_sum_ceiling']) / c['joint']['windows']
    for key, value in [('half_vector', c['vector']), ('charged_vector', c['joint']['vector'])]:
        if old[key] != value:
            raise ValueError('Fixed physical sources differ')
    if (F(old['ratio']), F(old['theta_half']), F(old['theta_joint']), F(old['W'])) != (ratio, th, tj, F(c['target']['W'])):
        raise ValueError('Fixed family differs')
    with gzip.open(args.ledger, 'rt') as handle:
        ledger = json.load(handle)
    count = len(ledger['states'])
    if not count == len(ledger['rows']) == len(ledger['energies']) <= 13000:
        raise ValueError('Bounded consistent candidate ledger required')
    rhs = [F(1)] + [F(0)] * 82 + [th, tj]
    if list(map(F, ledger['rhs'])) != rhs:
        raise ValueError('Ledger target differs')
    arr = np.array(ledger['rows']).T
    lp = linprog(ledger['energies'], A_eq=arr[:83], b_eq=np.array(rhs[:83], float),
                 A_ub=arr[83:], b_ub=np.array(rhs[83:], float), bounds=(0,None), method='highs',
                 options={'dual_feasibility_tolerance':1e-10,'primal_feasibility_tolerance':1e-10})
    if not lp.success:
        raise ValueError('No numerical candidate basis: ' + lp.message)
    values = list(lp.x) + list(np.array(rhs[83:],float) - arr[83:] @ lp.x)
    support = [i for i,v in enumerate(values) if v > 1e-14]
    n = len(support)
    if not 1<=85-n<=12:
        raise ValueError('Residual completion requires deficiency one through12')
    row_scale = lcm(hn, cn * ratio.denominator, th.denominator, tj.denominator)
    target = [1] + [0]*82 + [int(th*row_scale), int(tj*row_scale)]
    action = [actions()] + [spin_actions({k:1}) for k in LABELS] + [spectator_actions({k:1}) for k in SPECTATORS]

    def sparse(matrix):
        denominator = lcm(*(F(v).denominator for row in matrix for v in row))
        return denominator, [(i,j,int(F(v)*denominator)) for i,row in enumerate(matrix) for j,v in enumerate(row) if v]

    blocks = {}
    for key,a,ds,ph,q,ts,cols,k,derivatives in mats:
        blocks[key] = {
            'cols': cols, 'norms': [sum(v*v for v in col.values()) for col in cols],
            'ts': [[int(v) for v in row] for row in ts], 'energy': sparse(k),
            'coherent': [sparse(m) for m in derivatives] + [sparse(projected_matrix(cols,item)) for item in action],
        }
    cache = {}
    local=c['local_window'];onsite=list(map(F,local['onsite_profile']));density=list(map(F,local['density_profile']));range2=list(map(F,local['range_two_density_profile']))
    fixed_terms=[(quadratic_value,{k:F(v) for k,v in c['quadratic_charge_telescope'].items()}),(square_value,{k:F(v) for k,v in c['charge_square_pair_telescope'].items()}),(indicator_value,{k:F(v) for k,v in c['higher_charge_indicator_telescope'].items()})]

    def physical(vector_source):
        vector, norm = _physical_vector(vector_source)
        identity = tuple(sorted(vector.items()))
        if identity in cache:
            return cache[identity]
        if len(vector)==1:
            state,amplitude=next(iter(vector.items()))
            charges=[((state>>(2*i))&3).bit_count()-1 for i in range(6)]
            diagonal=[shape.get(state&1023,0)-shape.get(state>>2,0) for shape in shapes]+[int(signed_value(state,{label:F(1)})) for label in PATTERNS]
            ph=F(amplitude*amplitude*half.get(state,0)**2,hn)
            joint=ph+ratio*sum((F(amplitude*amplitude*source.get(state,0)**2,cn) for source in charged),F(0))
            row=[norm,0,0]+[norm*v for v in diagonal]+[int(norm*item[state].get(state,0)) for item in action]+[int(ph*row_scale),int(joint*row_scale)]
            energy=norm*(sum((a*q*q/2 for a,q in zip(onsite,charges)),F(0))+sum((a*charges[i]*charges[i+1] for i,a in enumerate(density)),F(0))+diagonal_value(state,range2)+sum((fn(state,terms) for fn,terms in fixed_terms),F(0)))
            result={'column':row,'energy':energy,'norm':norm,'vector':{str(state):amplitude}}
            cache[identity]=result
            return result
        sector = next(iter(vector))
        particles = ((sector & 1365).bit_count(), (sector & 2730).bit_count())
        found = None
        for key, block in blocks.items():
            if key[:2] != particles:
                continue
            coeff = [F(sum(a*vector.get(s,0) for s,a in col.items()),cnorm) for col,cnorm in zip(block['cols'],block['norms'])]
            if any(v.denominator != 1 for v in coeff):
                continue
            coeff = list(map(int, coeff))
            reconstructed = {s:z*a for z,col in zip(coeff,block['cols']) if z for s,a in col.items()}
            if reconstructed == vector:
                found = (block,coeff)
                break
        if found is None:
            raise ValueError('Physical vector outside reconstructed invariant columns')
        block, coeff = found
        if sum(z*z*k for z,k in zip(coeff,block['norms'])) != norm:
            raise ValueError('Physical norm mismatch')
        def expectation_raw(matrix):
            denominator, entries = matrix
            return F(sum(coeff[i]*coeff[j]*v for i,j,v in entries), denominator)
        coherent = [expectation_raw(matrix) for matrix in block['coherent']]
        if any(v.denominator != 1 for v in coherent):
            raise ValueError('Expected integral unit-operator moments')
        diagonal = [sum(v*z*z*k for v,z,k in zip(row,coeff,block['norms'])) for row in block['ts']]
        half_dot = sum(v*half.get(s,0) for s,v in vector.items())
        joint = F(half_dot*half_dot,hn) + ratio*sum((F(sum(v*source.get(s,0) for s,v in vector.items())**2,cn) for source in charged),F(0))
        row = [norm] + list(map(int,coherent[:2])) + diagonal + list(map(int,coherent[2:])) + [int(F(half_dot*half_dot,hn)*row_scale), int(joint*row_scale)]
        assert len(row) == 85
        energy = expectation_raw(block['energy']) - sum(p*d for p,d in zip(x[2:4],coherent[:2]))
        result = {'column':row,'energy':energy,'norm':norm,'vector':{str(s):v for s,v in vector.items()}}
        cache[identity] = result
        return result

    selected = []
    for index in support:
        if index >= count:
            selected.append({'column':[row_scale*int(j == 83+index-count) for j in range(85)],'energy':F(0),'norm':0,'vector':None})
        else:
            item = physical(ledger['states'][index]['vector'])
            normalized = np.array(item['column'],float)/item['norm']
            normalized[83:] /= row_scale
            if np.max(abs(normalized-arr[:,index])) > 1e-9:
                raise ValueError('Physical moments disagree with ledger')
            selected.append(item)
    additions = [physical(item['vector']) for item in old['mixture']]
    # Check the accepted source weights against the independently reconstructed
    # rows, with nonnegative fidelity slack, before using their columns.
    old_moments = [sum((F(item['weight'])*col['column'][i]/col['norm'] for item,col in zip(old['mixture'],additions)),F(0)) for i in range(85)]
    if old_moments[:83] != list(map(F,target[:83])) or any(v>b for v,b in zip(old_moments[83:],target[83:])):
        raise ValueError('Old source mixture failed independent moments')
    prepared = time.monotonic()
    d, transformed_rhs, transformed_old = reduce_basis([s['column'] for s in selected], target, [s['column'] for s in additions])
    reduced = time.monotonic()
    if not any(transformed_rhs[n:]):
        raise ValueError('Selected basis is already consistent; use the simpler replay')
    energy_scale = lcm(*(item['energy'].denominator for item in selected+additions))
    energies = [int(item['energy']*energy_scale) for item in selected+additions]
    # Solve the small residual system in physical-weight coordinates.
    # Normalize rhs and columns before numerical basis selection, then solve
    # that selected residual basis again with exact fractions.
    deficiency=85-n
    residual_rows=[[F(v,item['norm']) for v,item in zip(row,additions)] for row in transformed_old[n:]]
    row_scales=[max(map(abs,row)) for row in residual_rows]
    if min(row_scales)<=0:raise ValueError('Completion pool does not span residual rows')
    delta=max(F(abs(v))/scale for v,scale in zip(transformed_rhs[n:],row_scales))
    eq=np.array([[float(v/scale) for v in row] for row,scale in zip(residual_rows,row_scales)])
    eq_rhs=np.array([float(F(v)/scale/delta) for v,scale in zip(transformed_rhs[n:],row_scales)])
    inequalities=[];inequality_rhs=[]
    for i in range(n):
        row=[F(v,item['norm'])*delta for v,item in zip(transformed_old[i],additions)]
        bound=F(transformed_rhs[i])
        if d<0:row=[-v for v in row];bound=-bound
        scale=max(max(map(abs,row)),abs(bound))
        if scale:
            inequalities.append([float(v/scale) for v in row]);inequality_rhs.append(float(bound/scale))
    singular_values=np.linalg.svd(eq,compute_uv=False)
    rng=np.random.default_rng(20260912)
    def recover(indices):
        matrix=[[F(row[j]) for j in indices]+[F(value)] for row,value in zip(transformed_old[n:],transformed_rhs[n:])]
        pivot_row=0;pivots=[]
        for j in range(len(indices)):
            pivot=next((i for i in range(pivot_row,len(matrix)) if matrix[i][j]),None)
            if pivot is None:continue
            matrix[pivot_row],matrix[pivot]=matrix[pivot],matrix[pivot_row]
            scale=matrix[pivot_row][j];matrix[pivot_row]=[v/scale for v in matrix[pivot_row]]
            for i in range(len(matrix)):
                if i!=pivot_row and matrix[i][j]:
                    factor=matrix[i][j];matrix[i]=[v-factor*w for v,w in zip(matrix[i],matrix[pivot_row])]
            pivots.append(j);pivot_row+=1
        if any(row[-1] for row in matrix[pivot_row:]):raise ValueError('Dependent residual columns cannot represent exact target')
        free=[j for j in range(len(indices)) if j not in pivots]
        if len(free)>1:raise ValueError('More than one residual free variable requires another method')
        constants=[F(0)]*len(indices);slopes=[F(0)]*len(indices)
        for i,j in enumerate(pivots):
            constants[j]=matrix[i][-1]
            if free:slopes[j]=-matrix[i][free[0]]
        if not free:return constants,len(pivots)
        slopes[free[0]]=F(1)
        current_c=[(F(transformed_rhs[i])-sum((F(transformed_old[i][j])*v for j,v in zip(indices,constants)),F(0)))/d for i in range(n)]
        current_s=[-sum((F(transformed_old[i][j])*v for j,v in zip(indices,slopes)),F(0))/d for i in range(n)]
        lo=F(0);hi=None
        for constant,slope in zip(constants+current_c,slopes+current_s):
            if slope>0:lo=max(lo,-constant/slope)
            elif slope<0:hi=min(hi,-constant/slope) if hi is not None else -constant/slope
            elif constant<0:raise ValueError('No nonnegative affine completion')
        if hi is not None and hi<lo:raise ValueError('Empty exact nonnegative parameter interval')
        cost_slope=sum((v*item['energy'] for v,item in zip(current_s,selected)),F(0))+sum((v*additions[j]['energy'] for j,v in zip(indices,slopes)),F(0))
        if cost_slope<0 and hi is None:raise ValueError('Unexpected unbounded physical completion')
        parameter=hi if cost_slope<0 else lo
        return [c+s*parameter for c,s in zip(constants,slopes)],len(pivots)
    best=None;tried=feasible=bounded=0;attempts=[]
    for attempt in range(16):
        objective=np.ones(85) if attempt==0 else rng.uniform(.5,1.5,size=85)
        result=linprog(objective,A_eq=eq,b_eq=eq_rhs,A_ub=np.array(inequalities),b_ub=np.array(inequality_rhs),bounds=(0,None),method='highs',options={'dual_feasibility_tolerance':1e-10,'primal_feasibility_tolerance':1e-10})
        tried+=1;entry={'attempt':attempt,'status':result.message}
        if not result.success:attempts.append(entry);continue
        indices=[j for j,v in enumerate(result.x) if v>1e-12];entry['added_columns']=len(indices)
        if not 1<=len(indices)<=deficiency:entry['failure']='Residual basis exceeds missing dimension';attempts.append(entry);continue
        try:
            z,rank=recover(indices);entry['residual_rank']=rank
            if min(z)<0:raise ValueError('Negative exact added coefficient')
            current=[(F(transformed_rhs[i])-sum((F(transformed_old[i][j])*v for j,v in zip(indices,z)),F(0)))/d for i in range(n)]
            if min(current)<0:raise ValueError('Negative exact selected coefficient')
            feasible+=1;raw=current+z;columns=selected+[additions[j] for j in indices]
            weights=[v*item['norm'] for v,item in zip(raw,columns) if item['vector'] is not None]
            if any(len(str(v))>4096 or not 0<=v<=1 for v in weights):raise ValueError('Completed physical weight outside unchanged bounds')
            upper=sum((v*item['energy'] for v,item in zip(raw,columns)),F(0))/5
            entry.update(exact_feasible=True,upper_float=float(upper))
            if best is None or upper<F(best['cost_num'],best['cost_den']):
                bounded+=1;best={'indices':indices,'raw':raw,'columns':columns,'cost_num':upper.numerator,'cost_den':upper.denominator}
        except ValueError as error:entry['failure']=str(error)
        attempts.append(entry)
    diagnostic = {'accepted':False,'proposal_written':False,'selected_columns':n,'deficiency':85-n,
                  'completions_tried':tried,'nonnegative_completions':feasible,'improving_bounded_completions':bounded,
                  'preparation_seconds':prepared-started,'integer_reduction_seconds':reduced-prepared,
                  'pivot_bits':abs(d).bit_length(),'ledger_candidates':count,'completion_method':'bounded residual LP then exact solve','residual_dimension':deficiency,'residual_variable_scale':str(delta),'residual_singular_values':singular_values.tolist(),'residual_attempts':attempts}
    if best is not None:
        columns,raw = best['columns'],best['raw']
        if any(sum((v*item['column'][i] for v,item in zip(raw,columns)),F(0)) != target[i] for i in range(85)):
            raise ValueError('Completed basis failed original exact equations')
        mixture = [{'weight':str(v*item['norm']),'vector':item['vector']} for v,item in zip(raw,columns) if v and item['vector'] is not None]
        if len(mixture)>85 or sum(F(item['weight']) for item in mixture)!=1:
            raise ValueError('Physical mixture cap or trace failed')
        upper = F(best['cost_num'],best['cost_den'])
        proposal = {k:old[k] for k in ['half_vector','charged_vector','ratio','theta_half','theta_joint','diagonal_shapes','W','range_two_density_profile']}
        proposal.update(kind='joint_spectator_hopping_family_proposal_v1',mixture=mixture,proposed_periodic_family_upper=str(upper),scope='Untrusted exact completed-basis proposal; independent physical replay required.')
        (out/'diagonal_family_limit_proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
        diagnostic.update(proposal_written=True,mixture_sources=len(mixture),rational_family_upper=str(upper),rational_family_upper_float=float(upper),added_old_columns=best['indices'])
    files={Path(__file__).resolve(),Path(__file__).with_name('fraction_free_completion.py'),cp,old_path,old_receipt_path,args.ledger.resolve(),candidate_path}
    for module in tuple(sys.modules.values()):
        source=getattr(module,'__file__',None)
        if source and str(Path(source).resolve()).startswith(str(ROOT/'experiments')+'/'):
            files.add(Path(source).resolve())
    diagnostic['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
    diagnostic['seconds']=time.monotonic()-started
    (out/'basis_completion_diagnostic.json').write_text(json.dumps(diagnostic,indent=2)+'\n')
    print(json.dumps({k:v for k,v in diagnostic.items() if k not in ('source_sha256','rational_family_upper')}),flush=True)


if __name__ == '__main__':
    main()
