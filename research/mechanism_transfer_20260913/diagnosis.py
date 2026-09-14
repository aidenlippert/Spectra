"""Exact, non-overlapping relations behind the frozen H6 interval."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

ROOT=Path(__file__).resolve().parents[2]


def run():
    prior=ROOT/'results/trace_pricing_20260913';out=ROOT/'results/mechanism_transfer_20260913'
    audit=json.loads((prior/'audit.json').read_text());control=audit['physical_interval']
    receipt=json.loads((prior/'campaign/continue/round_10_batch_1/receipt.json').read_text())
    cert=json.loads((prior/'campaign/continue/round_10_batch_1/certificate.json').read_text())
    U=F(control['upper']);P=F(control['lower']);L=F(receipt['accepted']['original_lower_Ha'])
    C=F(audit['full_completed_spin_dual']['original_lower_ceiling_Ha'])
    tail=receipt['accepted']['tail'];ell=F(tail['lower_operator_shift_Ha']);u=F(tail['upper_operator_shift_Ha'])
    rho=F(receipt['accepted']['retained']['residual_l1']);b=F(cert['b'])
    if L!=b-rho+ell: raise AssertionError('Unexpected certificate decomposition')
    values={'interval_width':U-L,'upper_error_maximum':U-P,
        'actual_lower_error_minimum':P-L,'actual_lower_error_maximum':U-L,
        'family_true_error_floor':P-C,'fixed_upper_family_width_floor':U-C,
        'unresolved_distance_from_found_lower_to_ceiling':C-L,
        'tail_lower_shift_cost':-ell,'tail_combined_representation_envelope':u-ell,
        'exact_residual_penalty':rho,'retained_ground_minus_certificate_scalar_minimum':U-L-(U-P)-(u-ell)-rho,
        'export_loss_relative_to_floating_retained_objective':F(str(receipt['floating_lower_Ha']))-(L-ell)}
    result={'exact_Ha':{k:str(v) for k,v in values.items()},'mHa':{k:float(1000*v) for k,v in values.items()},
        'audit_sha256':hashlib.sha256((prior/'audit.json').read_bytes()).hexdigest(),
        'identity':'W=(U-E0(H))+(E0(H)-E0(Hret)-ell)+(E0(Hret)-b)+rho',
        'caveats':['The first two unknown terms lie respectively in [0,U-P] and [0,u-ell].',
            'b is the exact exported scalar before residual penalty, not an independently accepted lower.',
            'The export loss relative to the floating objective is a diagnostic comparison, not a second additive penalty.',
            'C need not be the optimal family ceiling. C-L cannot be assigned entirely to direction selection or entirely to family inadequacy.',
            'Representation, upper uncertainty and residual are bounded contributions; their simultaneous worst cases need not occur.']}
    (out/'interval_diagnosis.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':run()
