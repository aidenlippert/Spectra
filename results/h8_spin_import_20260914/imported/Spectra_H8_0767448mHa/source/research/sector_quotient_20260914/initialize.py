"""Freeze the request, inherited bytes, endpoint, and permitted dependencies."""
import json
from pathlib import Path
from research.reconstruction_compression_20260914.inputs import sha, dump
from research.sector_quotient_20260914.budget import ROOT, OUT


def run():
    OUT.mkdir(exist_ok=False)
    parent = ROOT / 'results/collective_completion_20260914'
    inherited = json.loads((parent/'preservation_before.json').read_text())['files']
    inherited.update(json.loads((parent/'manifest.json').read_text())['files'])
    for name in ('manifest.json', 'seal.json'):
        p = parent/name
        if p.exists():
            inherited[str(p.relative_to(ROOT))] = {'bytes': p.stat().st_size, 'sha256': sha(p)}
    errors = [p for p, rec in inherited.items() if not (ROOT/p).is_file() or sha(ROOT/p) != rec['sha256']]
    if errors:
        raise ValueError(('Changed inherited files', errors))
    dump(OUT/'preservation_before.json', {'files': inherited, 'checked_files': len(inherited)})
    frozen = json.loads((parent/'frozen_inputs.json').read_text())
    dump(OUT/'frozen_inputs.json', frozen)
    sources = [Path('/Users/aidenlippert/DERIVATION.md'), Path('/Users/aidenlippert/.codex/attachments/2a7df423-d724-49e6-9958-0d2a3e4099c9/pasted-text.txt')]
    dump(OUT/'request_sources.json', {'files': {str(p): sha(p) for p in sources},
        'request': str(sources[1]), 'supporting_derivation': str(sources[0]),
        'supplied_prototype_code_found': False,
        'parent_manifest_sha256': sha(parent/'manifest.json'),
        'strong_proof_use': 'diagnosis only; no strong factors in compact span selection',
        'compact_span_source': 'preserved paired H8 dual48 construction',
        'many_body_enumeration_authorized_for_discovery': False})
    print(json.dumps({'preserved_files': len(inherited), 'new_output': str(OUT)}))


if __name__ == '__main__':
    run()
