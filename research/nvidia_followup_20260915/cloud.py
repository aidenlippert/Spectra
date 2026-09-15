"""Own one explicitly authorized GPU experiment; never terminate other instances."""
import argparse
import datetime
import json
from pathlib import Path
from research.certificate_scaling.lambda_pool import _json

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/nvidia_followup_20260915'


def save(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT/name).open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def filtered(row):
    return {k: row[k] for k in ('id', 'name', 'status', 'ip', 'region', 'instance_type') if k in row}


def launch():
    if (OUT/'instance.json').exists():
        raise ValueError('Experiment already owns an instance')
    save('initial_instances.json', [filtered(r) for r in _json('GET', '/instances')['data']])
    kinds = _json('GET', '/instance-types')['data']
    typ = 'gpu_1x_a100_sxm4'
    spec = kinds[typ]
    price = spec['instance_type']['price_cents_per_hour']
    regions = {r['name'] for r in spec['regions_with_capacity_available']}
    if price > 199 or 'us-west-2' not in regions:
        raise ValueError('Recheck price or capacity before launching the intended A100')
    names = {row['name'] for row in _json('GET', '/ssh-keys')['data']}
    if 'aiden-mac' not in names:
        raise ValueError('Existing requested SSH key is unavailable')
    body = {'name': 'spectra-nvidia-h10-20260915', 'region_name': 'us-west-2',
        'instance_type_name': typ, 'ssh_key_names': ['aiden-mac'], 'quantity': 1}
    response = _json('POST', '/instance-operations/launch', body)['data']
    if len(response['instance_ids']) != 1:
        raise ValueError('Unexpected launch result; inspect provider state')
    record = {'id': response['instance_ids'][0], 'request': body, 'price_cents_per_hour': price,
        'launch_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat()}
    save('instance.json', record)
    print(json.dumps(record))


def status():
    owned = json.loads((OUT/'instance.json').read_text())['id']
    rows = _json('GET', '/instances')['data']
    row = next((filtered(r) for r in rows if r['id'] == owned), {'id': owned, 'status': 'absent'})
    if row.get('ip') and not (OUT/'ready.json').exists():
        save('ready.json', row)
    print(json.dumps(row))


def terminate():
    record = json.loads((OUT/'instance.json').read_text())
    rows = _json('GET', '/instances')['data']
    row = next((r for r in rows if r['id'] == record['id']), None)
    if row and row['name'] != record['request']['name']:
        raise ValueError('Owned-instance identity changed')
    response = _json('POST', '/instance-operations/terminate', {'instance_ids': [record['id']]}) if row else {'already_absent': True}
    save('termination.json', {'id': record['id'], 'requested_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'response': response, 'other_instances_untouched': True})
    print(json.dumps({'termination_requested_for_owned_instance': record['id']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('launch', 'status', 'terminate'))
    a = p.parse_args()
    {'launch': launch, 'status': status, 'terminate': terminate}[a.action]()
