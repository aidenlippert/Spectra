"""Exact joint/telescoping lower and matched physical upper on a transfer target."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_projector_extendibility import replay
from experiments.marginal_boundary_transfer import compile_block,contract,enclose
BASE=ROOT/'results/marginal_graded_hubbard8'


def main():
    parser=argparse.ArgumentParser();parser.add_argument('certificate',type=Path);args=parser.parse_args()
    started=time.monotonic();c=json.loads(args.certificate.read_text());target=c['target']
    source_path=BASE/'density_transfer/certificate.json';source=json.loads(source_path.read_text())
    choices=[item for item in source['targets'] if item['target']==target]
    if len(choices)!=1:raise ValueError('No matched physical transfer recipe for this target')
    case=choices[0]
    if c['chain_sites']!=source['sites'] or c['chain_sites']%8:raise ValueError('Matching multiple-of-eight chain required')
    lower=replay(c);print('exact lower accepted',flush=True)
    compiled=compile_block(source['upper'],source['hamiltonian'],**{k:F(v) for k,v in target.items()})
    exact=contract(compiled,case['a'],case['b'],3);small=enclose(compiled,case['a'],case['b'],3,160)
    if not F(small['energy_lower'])<=F(exact['energy'])<=F(small['energy_upper']):raise ValueError('Independent24-site contraction not enclosed')
    upper=enclose(compiled,case['a'],case['b'],c['chain_sites']//8,160)
    if upper['target']!=target or upper['sites']!=c['chain_sites']:raise ValueError('Mismatched physical upper')
    lo=F(lower['open_lower_density']);hi=F(upper['upper_per_site'])
    if lo>hi:raise ValueError('Lower exceeds physical upper')
    previous_path=BASE/'density_transfer/independent_replay.json';previous=json.loads(previous_path.read_text())
    prev=next(v for v in previous['targets'] if v['target']==target)
    files={Path(__file__).resolve(),args.certificate.resolve(),source_path,previous_path}
    for mod in tuple(sys.modules.values()):
        src=getattr(mod,'__file__',None)
        if src and str(Path(src).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(src).resolve())
    receipt={'accepted':True,'target':target,'sites':c['chain_sites'],'lower_per_site':str(lo),'upper_per_site':str(hi),
             'width_per_site':str(hi-lo),'previous_lower_per_site':prev['lower_per_site'],
             'lower_improvement':str(lo-F(prev['lower_per_site'])),'lower_replay':lower,'upper_replay':upper,
             'exact_24_site':exact,'enclosed_24_site':small,'seconds':time.monotonic()-started,
             'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
             'scope':'Fresh all-Fock joint/telescoping lower and physical transfer-state upper for the specified nearest-neighbor U,t,V chain. Projector sources and eight correction shapes came from V=+1/2; local coefficients were retuned for this target. Parameter transfer within one dimension, not a generic molecular or universal-representability result.'}
    args.certificate.with_name('matched_transfer_replay.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('lower_replay','upper_replay','exact_24_site','enclosed_24_site','source_sha256')}),flush=True)


if __name__=='__main__':main()
