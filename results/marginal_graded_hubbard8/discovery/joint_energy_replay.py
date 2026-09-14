"""Standard-library accepting replay for an explicit joint energy certificate."""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, sys, time
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_projector_extendibility import replay
BASE=ROOT/'results/marginal_graded_hubbard8'


def main():
    parser=argparse.ArgumentParser();parser.add_argument('certificate',type=Path);parser.add_argument('--psd-witnesses',type=Path)
    args=parser.parse_args();start=time.monotonic()
    certificate=json.loads(args.certificate.read_text())
    witnesses=json.loads(args.psd_witnesses.read_text())['witnesses'] if args.psd_witnesses is not None else None
    result=replay(certificate,psd_witnesses=witnesses)
    if result['local_sum_dimensions']!=4096 or len(result['local_sectors'])!=94:
        raise ValueError('Incomplete six-site all-Fock replay')
    previous_path=BASE/'six_site_projector/refined_independent_replay.json'
    previous=json.loads(previous_path.read_text())
    if result['chain_sites']!=previous['sites'] or result['target']!=previous['target']:
        raise ValueError('Historical comparison target mismatch')
    lower=F(result['open_lower_density']);upper=F(previous['upper_per_site'])
    if lower>upper:raise ValueError('New lower contradicts accepted physical upper')
    files={Path(__file__).resolve(),args.certificate.resolve(),previous_path}
    if args.psd_witnesses is not None:files.add(args.psd_witnesses.resolve())
    for module in tuple(sys.modules.values()):
        source=getattr(module,'__file__',None)
        if source and str(Path(source).resolve()).startswith(str(ROOT/'experiments')+'/'):
            files.add(Path(source).resolve())
    receipt={'accepted':True,'lower_replay':result,'lower_per_site':str(lower),
             'previous_lower_per_site':previous['lower_per_site'],
             'improvement_per_site':str(lower-F(previous['lower_per_site'])),
             'historical_upper_per_site':str(upper),'width_to_historical_upper':str(upper-lower),
             'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
             'seconds':time.monotonic()-start,
             'scope':'Fresh standard-library exact combined energy replay. Historical matched physical upper is referenced, not recomputed by this driver. Finite Hubbard construction only; no general chemistry or cost-versus-accuracy claim.'}
    suffix='_accelerated_replay.json' if args.psd_witnesses is not None else '_replay.json'
    out=args.certificate.with_name(args.certificate.stem+suffix)
    out.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('lower_replay','source_sha256')}),flush=True)


if __name__=='__main__':main()
