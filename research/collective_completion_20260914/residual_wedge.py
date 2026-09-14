"""Reuse the existing exact exterior-power bound on a checked residual."""
from fractions import Fraction as F
from research.certificate_scaling.wedge_spectral_bound import replay_residual


def improve(core,residual,receipt,witness):
    spectral=replay_residual(core,residual,witness)
    previous=F(receipt['lower']);lower=max(previous,F(spectral['lower']))
    receipt.update(coefficient_residual_lower=str(previous),residual_wedge=spectral,
                   lower=str(lower),lower_float=float(lower))
    return receipt


if __name__=='__main__':
    import argparse,hashlib,json,time
    from pathlib import Path
    from experiments.marginal_symbolic import verified_residual
    from research.certificate_scaling.spin_twirl import twirl
    from research.certificate_scaling.wedge_spectral_bound import propose
    p=argparse.ArgumentParser();p.add_argument('certificate');p.add_argument('out');a=p.parse_args()
    start=time.monotonic();path=Path(a.certificate);raw=path.read_bytes();cert=json.loads(raw)
    residual,_=verified_residual(cert['core'])
    if cert.get('spin_twirl'):residual=twirl(residual)
    cert['residual_wedge_witness']=propose(residual,cert['modes'],cert['particles'])
    out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
    (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
    receipt={'source_certificate':str(path),'source_sha256':hashlib.sha256(raw).hexdigest(),
             'postprocess_seconds':time.monotonic()-start,'many_body_states_enumerated':0,
             'status':'requires_independent_exact_replay'}
    (out/'postprocess.json').write_text(json.dumps(receipt,indent=2)+'\n')
    discovery=path.parent/'discovery.json'
    if discovery.exists():
        d=json.loads(discovery.read_text());d.update(residual_postprocess=receipt,
            certificate_bytes=(out/'certificate.json').stat().st_size,
            total_seconds=d['total_seconds']+receipt['postprocess_seconds'])
        (out/'discovery.json').write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)
