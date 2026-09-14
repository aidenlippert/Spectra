"""Untrusted scalar-threshold proposal at a new target; all corrections frozen."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import sys
import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_charged_projectors import charged_vectors
from experiments.marginal_range_two_density import diagonal_value
from experiments.marginal_quadratic_charge_telescope import local_value as quadratic
from experiments.marginal_charge_square_pairs import local_value as squares
from experiments.marginal_charge_indicator_telescope import local_value as indicators
from experiments.marginal_signed_charge_telescope import local_value as signed
from experiments.marginal_spin_word_telescope import local_value as words
from experiments.marginal_hopping_telescope import actions as hop_actions, projected_matrix
from experiments.marginal_spin_telescope import actions as spin_actions
from experiments.marginal_spectator_hopping import actions as spectator_actions
from experiments.marginal_pair_transfer import actions as pair_actions
from experiments.marginal_two_spectator_hopping import actions as two_actions
from experiments.marginal_three_spectator_hopping import actions as three_actions
from experiments.marginal_coherent_projector_telescope import actions as coherent_actions


def minimum(c):
    local = c['local_window']
    local_data = sector_matrices(6,F(local['U']),F(local['t']),
        list(map(F,local['onsite_profile'])),list(map(F,local['hopping_profile'])),
        F(local['V']),list(map(F,local['density_profile'])))
    action = [(F(c['hopping_telescope']),hop_actions())]
    for key,fn in [('spin_telescope',spin_actions),('spectator_hopping',spectator_actions),('pair_transfer',pair_actions),('two_spectator_hopping',two_actions),('three_spectator_hopping',three_actions),('coherent_projector',coherent_actions)]:
        if key in c:action.append((F(1),fn({k:F(v) for k,v in c[key].items()})))
    half = {int(s):a for s,a in c['vector'].items()}; hn = sum(a*a for a in half.values())
    charged,cn = charged_vectors(c['joint']['vector'])
    telescope = {int(s):F(a) for s,a in c['telescoping_diagonal'].items()}
    terms = [(fn,{k:F(v) for k,v in c[key].items()}) for key,fn in
             [('quadratic_charge_telescope',quadratic),('charge_square_pair_telescope',squares),
              ('higher_charge_indicator_telescope',indicators),('signed_charge_telescope',signed)]]
    terms.append((words,{k:F(v) for k,v in c.get('spin_word_telescope',{}).items()}))
    lowest = float('inf')
    for key,(_,matrix,cols) in local_data.items():
        scale = np.sqrt([sum(a*a for a in col.values()) for col in cols])
        denom = scale[:,None]*scale[None,:]
        a = np.array(matrix,float)/denom
        for weight,images in action:
            a += float(weight)*np.array(projected_matrix(cols,images),float)/denom
        diagonal = []
        for col in cols:
            values = {diagonal_value(s,list(map(F,local['range_two_density_profile']))) +
                      telescope.get(s&1023,F(0))-telescope.get(s>>2,F(0)) +
                      sum(fn(s,t) for fn,t in terms) for s in col}
            if len(values) != 1: raise ValueError('Lost reflection closure')
            diagonal.append(float(values.pop()))
        a += np.diag(diagonal)
        for vector,norm,penalty in [(half,hn,F(c['penalty'])+F(c['joint']['penalty']))] + [
            (v,cn,F(c['joint']['penalty'])*F(c['joint']['ratio'])) for v in charged]:
            w = np.array([sum(v*vector.get(s,0) for s,v in col.items()) for col in cols],float)/scale/np.sqrt(float(norm))
            a += float(penalty)*np.outer(w,w)
        lowest = min(lowest,float(eigh(a,subset_by_index=[0,0],eigvals_only=True)[0]))
    return lowest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seed',type=Path); parser.add_argument('output',type=Path)
    args = parser.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    c = json.loads(args.seed.read_text())
    original_minimum = minimum(c)
    if not 0 <= original_minimum-float(F(c['penalized_lower'])) < 1e-6:
        raise ValueError('Independent source spectrum does not match the rounded threshold')
    old_target = dict(c['target']); new_target = dict(U='5',t='1',V='1/4',W='-1/5')
    local = c['local_window']
    local['U'] = str(F(5,6)*F(new_target['U']))
    local['onsite_profile'] = [str(F(v)*F(new_target['U'])/F(old_target['U'])) for v in local['onsite_profile']]
    local['V'] = new_target['V']
    local['density_profile'] = [str(F(v)*F(new_target['V'])/F(old_target['V'])) for v in local['density_profile']]
    local['range_two_density_profile'] = [str(F(v)+F(5,4)*(F(new_target['W'])-F(old_target['W']))) for v in local['range_two_density_profile']]
    c['target'] = new_target
    value = minimum(c)
    c['penalized_lower'] = str(F(math.floor(value*10**7)-1,10**7))
    (args.output/'profile_joint_r1_2_certificate.json').write_text(json.dumps(c,indent=2)+'\n')
    lower = (F(c['penalized_lower'])-F(c['penalty'])*F(c['projector_sum_ceiling'])/c['windows']-
             F(c['joint']['penalty'])*F(c['joint']['projector_sum_ceiling'])/c['joint']['windows'])/5
    files = {Path(__file__).resolve(),args.seed.resolve()}
    for module in tuple(sys.modules.values()):
        path = getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'): files.add(Path(path).resolve())
    result = {'accepted':False,'old_target':old_target,'target':new_target,
              'source_minimum':original_minimum,'new_minimum':value,'proposed_periodic_lower':str(lower),
              'proposed_periodic_lower_float':float(lower),
              'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
              'scope':'New U,V,W target on the same chain geometry. All projector sources, penalties and correction coefficients frozen. Physical profiles rescaled to the new interactions; only the scalar local threshold recomputed. Untrusted proposal; no optimization, general chemistry or scalability claim.'}
    (args.output/'frozen_target_proposal.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'}))


if __name__ == '__main__': main()
