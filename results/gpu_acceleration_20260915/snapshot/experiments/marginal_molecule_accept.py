"""Compact the molecular coefficients and regenerate a tight Rayleigh witness."""
from fractions import Fraction as F
import json
import math
from pathlib import Path

from experiments.marginal_symbolic import decode, encode
from experiments.marginal_transfer_verify import replay


def run():
    # Numerical diagonalization only proposes amplitudes. replay recomputes the
    # energy exactly from the integer vector and exported Hamiltonian.
    from experiments.marginal_transfer import upper_ed
    root=Path(__file__).resolve().parents[1];out=root/'results/marginal_molecule'
    original=json.loads((out/'h4_rectangle_sto3g.json').read_text())
    cert=json.loads((out/'h4_degree3_certificate.json').read_text())
    old=decode(original['hamiltonian'],8,4)
    if old!=decode(cert['hamiltonian'],8,4):raise ValueError('Molecular fixture/certificate mismatch')
    denominator=10**12
    rounded={w:F(round(c*denominator),denominator) for w,c in old.items()}
    rounded={w:c for w,c in rounded.items() if c}
    extra=sum(abs(rounded.get(w,0)-c) for w,c in old.items())
    # Original generator stored float(exact_fraction_error). The next float
    # above it is a conservative rational enclosure of that rounded value's
    # preimage under nearest floating conversion.
    prior=F.from_float(math.nextafter(original['integral_rounding_l1_bound'],math.inf))
    budget=F(math.ceil((prior+extra)*10**15),10**15)
    fixture=dict(original);fixture['hamiltonian']=encode(rounded)
    fixture['additional_quantization_l1']=str(extra)
    fixture['numerical_integral_perturbation_l1_upper']=str(budget)
    fixture['coefficient_denominator']=denominator
    cert['hamiltonian']=encode(rounded)
    cert['independent_upper']=upper_ed(rounded,8,4)
    result=replay(cert)
    result['numerical_integral_perturbation_l1_upper']=str(budget)
    result['coefficient_perturbation_enclosed_lower']=str(F(result['lower'])-budget)
    result['coefficient_perturbation_enclosed_upper']=str(F(result['upper'])+budget)
    result['coefficient_denominator']=denominator
    result['scope']='Electronic energy of the rational finite-basis Hamiltonian; numerical integral and basis errors are not certified.'
    (out/'h4_compact_fixture.json').write_text(json.dumps(fixture,indent=2)+'\n')
    (out/'h4_compact_certificate.json').write_text(json.dumps(cert)+'\n')
    (out/'h4_compact_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({**result,'lower_float':float(F(result['lower'])),'upper_float':float(F(result['upper']))}),flush=True)
    return result


if __name__=='__main__':run()
