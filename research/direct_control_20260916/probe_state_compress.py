import json,time
from pathlib import Path
from .validated_tensor import *
s=json.loads(Path('results/direct_control_20260916/imported/Spectra_control_reduction/inputs/state.json').read_text())
a,e=load_mps(s);a,c=right_canonicalize(a)
for bond in [16,24,32,48,64,96,128]:
 t=time.monotonic();b,d,loc=compress(a,bond)
 print(json.dumps({'bond':bond,'error':e+c+d,'seconds':time.monotonic()-t,'norm':expectation(b)}),flush=True)
