import json,time
from pathlib import Path
from .validated_tensor import load_mps,zip_apply,expectation,tensor_bounds
p=Path('results/direct_control_20260916/imported/Spectra_control_reduction/inputs')
s=json.loads((p/'state.json').read_text());h=json.loads(Path('results/direct_control_20260916/exact_h8_mpo.json').read_text())
t=time.monotonic();a,e=load_mps(s);print('load',time.monotonic()-t,e,flush=True)
out,er,loc=zip_apply(a,h,64)
print(json.dumps({'seconds':time.monotonic()-t,'input_error':e,'action_error':er,'locals':loc,'bonds':[x.shape[2] for x in out]}),flush=True)
