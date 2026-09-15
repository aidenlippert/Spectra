import hashlib,json
from fractions import Fraction as F
from experiments.marginal_molecular_gershgorin import matrix

P='results/certificate_scaling/active_space_ladder/h4/fixture.json'
d=json.load(open(P)); states,a=matrix(d)
diag=[a[i][i] for i in range(70)]
rows=[diag[i]-sum(abs(a[i][j]) for j in range(70) if j!=i) for i in range(70)]
print({'fixture_sha256':hashlib.sha256(open(P,'rb').read()).hexdigest(),
       'diag_min':str(min(diag)),'gershgorin_min':str(min(rows)),
       'diag_min_float':float(min(diag)),'gershgorin_min_float':float(min(rows)),
       'derivation':'H_ii - sum_{j!=i}|H_ij| <= lambda_min(H)'} )
assert hashlib.sha256(open(P,'rb').read()).hexdigest()=='8e64505deb4cc06a87f75e669b0f3073d15ebca72daa6d086e06eeea0b0f5120'
assert min(rows) <= min(diag)
