"""Separate LP proposal error from rational rounding without weakening replay."""
from pathlib import Path
from fractions import Fraction as F
from unittest.mock import patch
import json,time
import experiments.marginal_joint_coefficient_constructor as joint
from experiments.marginal_polynomial_metric import JointPolynomial,check_positive

root=Path('results/marginal_graded_hubbard8');data=json.loads((root/'hamiltonian.json').read_text());proposal=json.loads((root/'joint_reynolds_conditional_extended/proposal.json').read_text());out=root/'conditional_export_diagnostic';out.mkdir(exist_ok=True)
original=joint.replay;records=[]
for name,precision,denominator in [('default',10**9,10**12),('fine',10**12,10**14)]:
    captured=[];started=time.monotonic()
    def audit(certificate):captured.append(certificate);return original(certificate)
    try:
        with patch.object(joint,'replay',audit):receipt=joint.export(data,proposal,out/name,metric_precision=precision,proof_denominator=denominator)
        record={'name':name,'accepted':True,'receipt':receipt}
    except ValueError as error:
        if not captured:raise
        certificate=captured[-1];ring=JointPolynomial(certificate);w,k,cost=ring.compile(F(certificate['target_lower']))
        parts={p:check_positive(ring,poly,scale,certificate[p+'_proof']) for p,poly,scale in [('weight',w,cost['metric_scale']),('numerator',k,cost['numerator_scale'])]}
        (out/(name+'_rejected_candidate.json')).write_text(json.dumps(certificate,indent=2)+'\n')
        record={'name':name,'accepted':False,'error':str(error),'parts':parts}
    record.update(metric_precision=precision,proof_denominator=denominator,seconds=time.monotonic()-started);records.append(record);print(json.dumps(record),flush=True)
(out/'diagnostic.json').write_text(json.dumps({'records':records,'scope':'Same numerical proposal at two rational precisions; full unchanged replay is the acceptance gate. Rejected candidates are retained only as diagnostics.'},indent=2)+'\n')
