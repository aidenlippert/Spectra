import json,time,argparse
from pathlib import Path
from .validated_tensor import *
p=argparse.ArgumentParser();p.add_argument('--bond',type=int,default=64);args=p.parse_args()
base=Path('results/direct_control_20260916/imported/Spectra_control_reduction/inputs')
s=json.loads((base/'state.json').read_text());h=json.loads(Path('results/direct_control_20260916/exact_h8_mpo.json').read_text())
t=time.monotonic();a,e=load_mps(s);a,ca=right_canonicalize(a)
w,we=dense_mpo(h);w,wc=right_canonicalize(w)
print(json.dumps({'preparation_seconds':time.monotonic()-t,'input_conversion_error':e,'state_canonicalization_error':ca,'operator_conversion_error':we,'operator_canonicalization_error':wc}),flush=True)
out,er,loc=zip_apply_dense(a,w,args.bond)
print(json.dumps({'seconds':time.monotonic()-t,'action_compression_error':er,'locals':loc,'bonds':[x.shape[2] for x in out]}),flush=True)
