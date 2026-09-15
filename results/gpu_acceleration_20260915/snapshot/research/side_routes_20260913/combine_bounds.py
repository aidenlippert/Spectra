"""Accept bounds from distinct positive amplitudes only for the same Hamiltonian."""
from fractions import Fraction as F
import hashlib
import json

from research.side_routes_20260913.finite_range import replay


def combine(payloads):
    if not payloads:
        raise ValueError('No candidate certificates')
    model = payloads[0]['model']
    if any(p['model'] != model for p in payloads):
        raise ValueError('Cannot combine bounds for different declared Hamiltonians')
    low_index = max(range(len(payloads)),key=lambda i:F(payloads[i]['claim']['lower']))
    high_index = min(range(len(payloads)),key=lambda i:F(payloads[i]['claim']['upper']))
    selected = {i:replay(payloads[i]) for i in {low_index,high_index}}
    lower, upper = F(selected[low_index]['lower']), F(selected[high_index]['upper'])
    if lower > upper:
        raise ValueError('Combined bounds are inconsistent')
    return {'lower':str(lower),'upper':str(upper),'width':str(upper-lower),
            'lower_float':float(lower),'upper_float':float(upper),'width_float':float(upper-lower),
            'lower_source_index':low_index,'upper_source_index':high_index,
            'selected_certificates_replayed':len(selected),
            'model_sha256':hashlib.sha256(json.dumps(model,sort_keys=True).encode()).hexdigest(),
            'selected_verification_seconds':sum(r['wall_seconds'] for r in selected.values()),
            'scope':'Exact combination of independently replayed lower and upper certificates for an identical declared model.'}
