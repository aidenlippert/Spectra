"""Frozen-rule fresh geometry generation with separated reference work."""
from fractions import Fraction as F
import hashlib
import json
import time

from research.global_response_20260913 import global_program as gp


def run():
    frozen=gp.OUT/'frozen_rule.json'
    rule=json.loads(frozen.read_text())
    source=gp.ROOT/'research/global_response_20260913/global_program.py'
    if hashlib.sha256(source.read_bytes()).hexdigest()!=rule['global_program_sha256']:
        raise ValueError('Selected program changed after the rule freeze')
    folder=gp.OUT/'fresh_h6_1p73'
    if folder.exists():raise ValueError('Preserve the existing fresh fixture')
    from research.mechanism_transfer_20260913.fixtures import build
    from research.molecular_collective_20260913.core import extract,propose_tail
    start=time.monotonic();generation=build(6,F(173,100),folder,fci_control=True)
    data=json.loads((folder/'fixture.json').read_text());p=extract(data)
    tail,stats=propose_tail(data,p,10);(folder/'tail.json').write_text(json.dumps(tail,indent=2)+'\n')
    (folder/'reference_upper.json').write_bytes((folder/'upper.json').read_bytes())
    receipt={'generation':generation,'tail_discovery':stats,'wall_seconds':time.monotonic()-start,
        'frozen_rule_sha256':hashlib.sha256(frozen.read_bytes()).hexdigest(),
        'FCI_determinants_enumerated_for_reference':generation['fci_determinant_dimension'],
        'reference_is_lower_discovery_input':'Only its scalar target; amplitudes are excluded from response/factor generation.',
        'scope':'New 1.73 Angstrom model after freezing the rule. The FCI reference is explicitly charged and keeps the full interval reference-assisted.'}
    (folder/'construction.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'seconds':receipt['wall_seconds'],'reference_dimension':generation['fci_determinant_dimension']}),flush=True)


if __name__=='__main__':run()
