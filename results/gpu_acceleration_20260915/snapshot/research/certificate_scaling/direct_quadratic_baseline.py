"""Full quadratic Gram baseline, H-only discovery for a fair cone comparison."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.marginal_coefficient import dictionaries, solve_coefficients, export
from experiments.marginal_symbolic import decode

if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--fixture',type=Path,required=True)
    p.add_argument('--outputdir',type=Path,required=True)
    a=p.parse_args()
    f=json.loads(a.fixture.read_text())
    h=decode(f['hamiltonian'],f['modes'],4)
    blocks=dictionaries(f['modes'],'quadratic')
    sol=solve_coefficients(h,f['modes'],f['particles'],blocks)
    cert,rec=export(h,f['modes'],f['particles'],blocks,sol)
    rec['source_factors_used']=False
    rec['source_upper_used']=False
    rec['full_Fock_space_constructed']=False
    rec['certificate_bytes']=len(json.dumps(cert,separators=(',',':')).encode())
    a.outputdir.mkdir(parents=True,exist_ok=True)
    (a.outputdir/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
    (a.outputdir/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps(rec))
