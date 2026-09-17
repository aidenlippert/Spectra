"""A fixed-identity residual-only limit, not a limitation of the Gram family."""
from fractions import Fraction
import hashlib
import json
from research.acceptance_channels_20260915.campaign import OUT,dump


def main():
    folder=OUT/'replays/dense_t2_long'
    paths={name:folder/name for name in ('original_interval.json','exact/lower.json','exact/certificate.json')}
    interval=json.loads(paths['original_interval.json'].read_text())
    lower=json.loads(paths['exact/lower.json'].read_text())
    width=Fraction(interval['upper_Ha'])-Fraction(interval['lower_Ha'])
    if width!=Fraction(interval['width_Ha']):raise ValueError('Inconsistent interval')
    eta=Fraction(lower['singlet']['residual_l1'])
    best_residual_only=width-2*eta
    result={'scope':'Keep this exported b, positive squares, upper and other allowances fixed. Replace only the scalar residual lower -eta by another scalar lower r.',
        'argument':'The existing L1 certificate gives -eta I <= R <= eta I in the singlet sector. Any valid scalar lower r satisfies r <= eta. Replacing -eta by r can therefore improve the lower by at most 2 eta. This does not bound improvements obtained by changing b, reoptimizing the squares, or using their joint spectrum with R.',
        'width_Ha':str(width),'residual_eta_Ha':str(eta),
        'residual_share_percent':float(100*eta/width),
        'optimistic_width_after_residual_only_replacement_Ha':str(best_residual_only),
        'optimistic_width_after_residual_only_replacement_mHa':float(1000*best_residual_only),
        'target_Ha':'1/625','residual_only_route_cannot_meet_target':best_residual_only>Fraction(1,625),
        'exact_family_obstruction':False,
        'source_sha256':{name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in paths.items()}}
    dump(OUT/'fixed_identity_residual_diagnosis.json',result)
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
