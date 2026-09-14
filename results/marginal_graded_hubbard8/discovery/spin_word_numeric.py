"""Numerical full diagonal basis; exact membership of the frozen legacy span."""
from fractions import Fraction as F
import numpy as np
from signed_charge_numeric import prepare as previous_prepare,BASE
from experiments.marginal_spin_word_telescope import PATTERNS,five_site_value,local_value
from experiments.marginal_signed_charge_telescope import five_site_value as charge_value

def legacy_terms(c):
    sparse={int(s):F(v) for s,v in c['telescoping_diagonal'].items()}
    signed={key:F(v) for key,v in c['signed_charge_telescope'].items()}
    values={s:sparse.get(s,F(0))+charge_value(s,signed) for s in range(1024)}
    terms={key:values[int(key)] for key in PATTERNS}
    if any(five_site_value(s,terms)!=values[s] for s in range(1024)):
        raise ValueError('Selected legacy diagonal is outside the full spin-word space')
    return terms

def prepare(c,shapes):
    x,physical,mats,*rest=previous_prepare(c,shapes)
    for shape in shapes:
        terms={key:shape.get(int(key),F(0)) for key in PATTERNS}
        if any(five_site_value(s,terms)!=shape.get(s,0) for s in range(1024)):
            raise ValueError('Selected legacy shape is outside the spin-word space')
    converted=[]
    for key,a,ds,ph,q,old_ts,cols,k,derivatives in mats:
        ts=np.array([[int(local_value(next(iter(col)),{label:F(1)})) for col in cols] for label in PATTERNS])
        converted.append((key,a,ds,ph,q,ts,cols,k,derivatives))
    return x,physical,converted,*rest
