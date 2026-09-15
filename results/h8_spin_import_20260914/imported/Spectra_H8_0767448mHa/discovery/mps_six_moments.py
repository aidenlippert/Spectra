from common import *
# Load the existing pure dagger helper without importing its unused CVXPY CLI.
import ast,types
source=ROOT/'experiments/marginal_coefficient.py'
node=next(n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='dagger')
shim=types.ModuleType('experiments.marginal_coefficient')
exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),shim.__dict__)
sys.modules['experiments.marginal_coefficient']=shim
from research.reconstruction_compression_20260914.moments import Proposal
from experiments.marginal_coefficient import dagger
start=time.monotonic();meta=json.loads((PREP/'frame.json').read_text())
fixture=json.loads((ROOT/'results/certificate_scaling/active_space_ladder/h8/fixture.json').read_text())
state=json.loads((ROOT/'results/correlated_pair_20260913/mps/h8_spatial_warm144/state.json').read_text())
words=[tuple(map(tuple,w)) for w in meta['rows']]
p=Proposal(fixture,state);print('START',len(words),'norm',p.norm,flush=True)
p.fill(words);print('FILL',time.monotonic()-start,'steps',p.steps,flush=True)
v=np.array([p.cache[w] for w in words]);np.save(OUT/'mps_six_moments.npy',v)
# Coefficient rows represent real-Hermitian polynomials; the stored representative
# off diagonal is paired with its adjoint in the map, so the coordinate factor is 2.
mul=np.array([1 if [i for c,i in w if c] == [i for c,i in w if not c] else 2 for w in words])
np.save(OUT/'mps_six_dual_coordinates.npy',v*mul)
record_json(OUT/'mps_six_moments.json',{'seconds':time.monotonic()-start,'words':len(words),'steps':p.steps,'norm':p.norm,'status':'proposal moments only, no energy acceptance','many_body_states_enumerated':0})
