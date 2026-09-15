"""Exact two-sided replay for frozen-H cubic certificate campaign artifacts."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import verify
from research.certificate_scaling.streaming_reference_upper import upper
from research.certificate_scaling.wedge_residual_bound import replay as wedge_replay
from research.certificate_scaling.wedge_spectral_bound import extract,replay_residual

def run(certificate,reference,out,method='coefficient',proof=None):
    start=time.monotonic();raw=certificate.read_bytes();cert=json.loads(raw);refraw=reference.read_bytes();ref=json.loads(refraw)
    sha=hashlib.sha256(raw).hexdigest()
    if method=='coefficient':lower_receipt=verify(cert)
    elif method=='wedge':lower_receipt=wedge_replay(cert)
    elif method=='spectral':
        if proof is None:raise ValueError('Spectral factors required')
        witness=json.loads(proof.read_text())
        if witness['certificate_sha256']!=sha:raise ValueError('Spectral witness source mismatch')
        residual,_=extract(cert);lower_receipt=replay_residual(cert,residual,witness)
    else:raise ValueError('Unknown exact lower method')
    lo=F(lower_receipt['lower']);hi,upper_receipt=upper(cert,ref['independent_upper']);width=hi-lo
    if width<0:raise AssertionError('Inconsistent independently replayed interval')
    rec={'method':method,'lower':str(lo),'upper':str(hi),'width':str(width),'lower_float':float(lo),'upper_float':float(hi),
         'width_float':float(width),'passes_0_0016_Ha':width<=F(16,10000),'modes':cert['modes'],'particles':cert['particles'],
         'certificate':str(certificate),'certificate_bytes':len(raw),'certificate_sha256':sha,'reference':str(reference),
         'reference_sha256':hashlib.sha256(refraw).hexdigest(),'proof':str(proof) if proof else None,
         'proof_sha256':hashlib.sha256(proof.read_bytes()).hexdigest() if proof else None,
         'proof_bytes':proof.stat().st_size if proof else 0,'lower_replay':lower_receipt,'upper_replay':upper_receipt,
         'wall_seconds':time.monotonic()-start,'scope':'Frozen rational Hamiltonian in fixed total-N sector. Reference FCI discovery remains an exponential validation cost; no basis/model/experimental error guarantee.'}
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps({k:v for k,v in rec.items() if k not in ['lower_replay','upper_replay']}),flush=True)
    return rec

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,required=True);p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--method',choices=['coefficient','wedge','spectral'],default='coefficient');p.add_argument('--proof',type=Path)
    a=p.parse_args();run(a.certificate,a.reference,a.out,a.method,a.proof)
