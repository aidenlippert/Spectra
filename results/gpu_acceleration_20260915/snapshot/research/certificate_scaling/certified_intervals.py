"""Attach existing physical upper witnesses and replay compressed intervals.

The lower-certificate compression never consults the witness. Upper replay
enumerates 70 fixed-N states: this validation cost is explicitly accounted for.
"""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.marginal_transfer_verify import replay


def main():
    out=ROOT/'results/certificate_scaling/intervals'
    out.mkdir(exist_ok=True)
    witness=json.loads((ROOT/'results/marginal_molecule_stress/accepted_degree3_certificate.json').read_text())['independent_upper']
    setups=[('square',ROOT/'results/lambda_runs/certificate_scaling/retrieved/results/certificate_scaling/lambda_batch/sparse_factors',10**9),
            ('rectangle',ROOT/'results/certificate_scaling/rectangle_sparse_transfer',10**7)]
    records=[]
    for geometry,folder,original_denominator in setups:
        for denominator in (original_denominator,10**5,10**6):
            source=folder/f'd{denominator}_t0-1.json'
            certificate=json.loads(source.read_text())
            if geometry=='square':
                certificate['independent_upper']=copy.deepcopy(witness)
            target=out/f'{geometry}_d{denominator}_interval_certificate.json'
            target.write_text(json.dumps(certificate,separators=(',',':'))+'\n')
            receipt=replay(certificate)
            records.append({'geometry':geometry,'denominator':denominator,
                'certificate':str(target.relative_to(ROOT)),'certificate_bytes':target.stat().st_size,
                'interval':receipt,'upper_witness_sector_states':70,
                'upper_witness_source':'existing integer amplitudes; Rayleigh value recomputed for this Hamiltonian',
                'target_hartree':'3/2000','passes_target':F(receipt['width'])<F(3,2000),
                'scope':'rational finite-basis electronic Hamiltonian; excludes continuum/integral-algorithm errors'})
    (out/'summary.json').write_text(json.dumps(records,indent=2)+'\n')
    for r in records:
        print(r['geometry'],r['denominator'],r['interval']['width_float'],r['passes_target'])


if __name__=='__main__':
    main()
