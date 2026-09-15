"""Independently replay both freshly constructed endpoints and the freeze."""
from fractions import Fraction as F
import json
from pathlib import Path
from research.correlated_pair_20260913.fresh_transfer import hashes
from research.correlated_pair_20260913.mps_exact import check as upper_check
from research.correlated_pair_20260913.self_consistent.fixed_guide import check as lower_check


def run():
    base=Path('results/correlated_pair_20260913');folder=base/'fresh_h6_1p91';freeze=json.loads((base/'transfer_freeze.json').read_text())
    if freeze['source_sha256']!=hashes():raise ValueError('Frozen construction rule changed')
    data=json.loads((folder/'fixture.json').read_text());state=json.loads((folder/'mps/state.json').read_text());sos=json.loads((folder/'structured/certificate.json').read_text())
    uc=upper_check(data,state);lc=lower_check(data,sos);U=F(uc['upper_Ha']);L=F(lc['lower_Ha'])
    if U<L:raise ValueError('Inconsistent transfer endpoints')
    result={'case':'new H6 at1.91 Angstrom','upper':uc,'lower':lc,'lower_Ha':str(L),'upper_Ha':str(U),
            'width_mHa':float(1000*(U-L)),'target_1p6mHa_met':U-L<=F(16,10000),
            'frozen_sources_unchanged':True,'FCI_or_full_cubic_teacher_built':False,
            'scope':'A transferred full interval from independent Hamiltonian-only construction; width is reported without implying a transferred tight interval.'}
    (folder/'interval.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('upper','lower')},indent=2))


if __name__=='__main__':run()
