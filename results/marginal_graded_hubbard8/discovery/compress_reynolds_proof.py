"""Exact Shannon merging of complementary occupation indicators.

I_(R,O) + I_(R,O xor bit) = I_(R without bit,O without bit).
Both nonnegative atom families admit this identity independently.
"""
from pathlib import Path
from collections import defaultdict
import copy,json,hashlib,time
from experiments.marginal_polynomial_metric import replay


def compress(atoms):
    counts=defaultdict(int)
    for atom in atoms:counts[atom['required'],atom['occupied']]+=atom['weight']
    merges=0
    for degree in range(max((r.bit_count() for r,o in counts),default=0),0,-1):
        for required in sorted({r for r,o in counts if r.bit_count()==degree}):
            for bit in (1<<i for i in range(required.bit_length()) if required&(1<<i)):
                for occupied in sorted(o for r,o in counts if r==required and not o&bit):
                    a=(required,occupied);b=(required,occupied|bit)
                    common=min(counts.get(a,0),counts.get(b,0))
                    if not common:continue
                    counts[a]-=common;counts[b]-=common;counts[required^bit,occupied]+=common;merges+=1
    result=[{'required':r,'occupied':o,'weight':weight} for (r,o),weight in sorted(counts.items()) if weight]
    return result,merges

if __name__=='__main__':
    root=Path('results/marginal_graded_hubbard8');source=root/'fixed_reynolds/certificate.json';raw=source.read_bytes();certificate=json.loads(raw);changes={};start=time.monotonic()
    for part in ('weight','numerator'):
        for family in ('positive_indicators','charge_indicators'):
            atoms=certificate[part+'_proof'][family];replacement,merges=compress(atoms)
            certificate[part+'_proof'][family]=replacement;changes[part+'/'+family]={'before':len(atoms),'after':len(replacement),'merges':merges}
    receipt=replay(certificate);out=root/'reynolds_compressed';out.mkdir(exist_ok=True)
    for name,value in [('certificate',certificate),('receipt',receipt),('compression',{'source_sha256':hashlib.sha256(raw).hexdigest(),'changes':changes,'seconds':time.monotonic()-start,'scope':'Exact indicator partition identities only; unchanged Hamiltonian, metric, bounds and number identities. Full original replay required.'})]:
        (out/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')
    print(json.dumps(changes),flush=True)
