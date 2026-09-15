"""Compressed overlapping three-mode CAR routing, in abstract energy units.

K_pq = n_p+n_q-2 n_p n_q - (a_p^dagger parity(p,q) a_q + h.c.).
For q=p+1 or p+2 this is a local occupation-swap Laplacian. We certify
sum J_i K_i,i+1 - sum t_i K_i,i+2 + a diagonal nearest-neighbor potential.
All verification and finite-grammar discovery arithmetic is rational.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import time

from research.side_routes_20260913.positive_chain import bounds as diagonal_bounds
from research.interference_routing_20260913.routing import rational

RATIOS = (F(1), F(1, 2), F(2), F(1, 4), F(4))


def validate(model):
    m, n = model['modes'], model['particles']
    if type(m) is not int or type(n) is not int or not 3 <= m <= 4096 or not 0 <= n <= m:
        raise ValueError('Invalid template sector/work budget')
    arrays = []
    for key, length in (('nearest', m-1), ('next_nearest', m-2), ('fields', m), ('interaction', m-1)):
        raw = model[key]
        if not isinstance(raw, list) or len(raw) != length:
            raise ValueError('Wrong array length: '+key)
        values = [rational(x) for x in raw]
        if key in ('nearest', 'next_nearest') and any(x < 0 for x in values):
            raise ValueError('Nonnegative routing magnitudes required')
        arrays.append(values)
    return m, n, *arrays, rational(model.get('offset', 0))


def discover(model):
    """Finite DP over five exact route ratios; failure is not impossibility."""
    start = time.monotonic()
    m, n, nearest, next_nearest, *_ = validate(model)
    previous = {F(0): None}
    layers = []
    transitions = 0
    peak = 1
    for i, t in enumerate(next_nearest):
        layer = {}
        for rho in RATIOS:
            a, b = t*(1+rho), t*(1+1/rho)
            for old_rho in previous:
                transitions += 1
                incoming = F(0) if i == 0 else next_nearest[i-1]*(1+1/old_rho)
                if incoming+a <= nearest[i] and (i < m-3 or b <= nearest[-1]):
                    layer[rho] = old_rho
                    break
        if not layer:
            raise ValueError(f'No route in five-ratio grammar at triangle {i}; not a general impossibility proof')
        layers.append(layer)
        previous = layer
        peak = max(peak, len(layer))
    rho = next(iter(previous))
    ratios = [None]*(m-2)
    for i in range(m-3, -1, -1):
        ratios[i] = str(rho)
        rho = layers[i][rho]
    return ratios, {'discovery_transitions': transitions, 'discovery_peak_states': peak,
                    'discovery_seconds': time.monotonic()-start, 'ratio_vocabulary': list(map(str, RATIOS))}


def replay(payload):
    start = time.monotonic()
    model = payload['model']
    m, n, nearest, next_nearest, fields, interaction, offset = validate(model)
    raw = payload['ratios']
    if not isinstance(raw, list) or len(raw) != m-2:
        raise ValueError('One ratio per triangle required')
    ratios = [rational(x) for x in raw]
    if any(x <= 0 for x in ratios):
        raise ValueError('Positive route ratios required')
    used = [F(0)]*(m-1)
    for i, (t, rho) in enumerate(zip(next_nearest, ratios)):
        # a K12+b K23-t K13 has capacity [[a-t,-t],[-t,b-t]].
        a, b = t*(1+rho), t*(1+1/rho)
        if a-t < 0 or b-t < 0 or (a-t)*(b-t) < t*t:
            raise ValueError('Local joint capacity is not PSD')
        used[i] += a
        used[i+1] += b
    if any(cost > available for cost, available in zip(used, nearest)):
        raise ValueError('Shared positive edge capacity overdrawn')
    diagonal = {'modes': m, 'particles': n, 'hopping': ['0']*(m-1),
                'interaction': list(map(str, interaction)), 'fields': list(map(str, fields)), 'offset': str(offset)}
    amplitude = {'sites': ['1']*m, 'bonds': ['1']*(m-1)}
    exact = diagonal_bounds(diagonal, amplitude)
    claims = {k: exact[k] for k in ('lower', 'upper', 'width')}
    if 'claim' in payload and payload['claim'] != claims:
        raise ValueError('False stored interval')
    return {**claims, 'lower_float': exact['lower_float'], 'upper_float': exact['upper_float'],
            'width_float': exact['width_float'], 'modes': m, 'particles': n,
            'routing_templates': m-2, 'largest_capacity_dimension': 2,
            'nonstoquastic_triangles': sum(t > 0 for t in next_nearest),
            'unused_positive_capacity': list(map(str, (a-b for a, b in zip(nearest, used)))),
            'diagonal_dp_transitions': exact['dp_transitions'], 'diagonal_peak_states': exact['peak_dp_states'],
            'peak_diagonal_dp_rational_bits': exact['peak_rational_bits'],
            'many_body_states_enumerated': 0, 'replay_seconds': time.monotonic()-start,
            'energy_units': 'abstract model unit, not Hartree',
            'scope': 'Declared density-assisted exchange family; uniform fixed-N upper; accuracy not generic.'}


def campaign(out):
    campaign_start = time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for m in (4, 8, 16, 32, 64):
        for name, delta in (('uniform', F(0)), ('alternating', F(0)),
                            ('weak_density', F(1, 10000)), ('strong_density', F(1, 5)),
                            ('varying_fields', F(0))):
            alternate = name == 'alternating'
            model = {'modes': m, 'particles': m//2,
                     'nearest': [str(F(3, 10) if i%2 == 0 else F(3, 5)) if alternate else '1' for i in range(m-1)],
                     'next_nearest': ['1/10' if alternate else '1/4']*(m-2),
                     'fields': [str(F((i%3)-1, 10)) if name == 'varying_fields' else '0' for i in range(m)],
                     'interaction': [str(delta)]*(m-1), 'offset': '0'}
            ratios, cost = discover(model)
            payload = {'model': model, 'ratios': ratios}
            result = replay(payload)
            payload['claim'] = {k: result[k] for k in ('lower', 'upper', 'width')}
            filename = f'{name}_{m}.json'
            (out/filename).write_text(json.dumps(payload, separators=(',', ':'))+'\n')
            ablation_start = time.monotonic()
            uniform_status = 'accepted'
            try:
                replay({'model': model, 'ratios': ['1']*(m-2)})
            except ValueError as exc:
                uniform_status = str(exc)
            rows.append({'case': name, 'file': filename, 'receipt': result, 'discovery': cost,
                         'all_equal_route_ablation': uniform_status,
                         'ablation_seconds': time.monotonic()-ablation_start,
                         'witness_bytes': (out/filename).stat().st_size})
    summary = {'cases': rows, 'campaign_seconds_before_summary_write': time.monotonic()-campaign_start,
               'accounting_scope': 'All 25 discoveries, accepting replays, equal-ratio ablations, proof serialization and file writes. Excludes interpreter startup and final summary serialization/write.'}
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({'cases': len(rows), 'largest_modes': 64, 'global_states_enumerated': 0}))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--campaign', type=Path)
    p.add_argument('--proof', type=Path)
    p.add_argument('--out', type=Path)
    args = p.parse_args()
    if args.campaign:
        campaign(args.campaign)
    else:
        receipt = replay(json.loads(args.proof.read_text()))
        args.out.write_text(json.dumps(receipt, indent=2)+'\n')
        print(json.dumps(receipt))
