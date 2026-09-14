from fractions import Fraction as F
import json, os
import numpy as np

def eig2(a):
    # exact characteristic polynomial for symmetric 2x2
    tr=a[0][0]+a[1][1]; det=a[0][0]*a[1][1]-a[0][1]*a[1][0]
    disc=tr*tr-4*det
    return (float(tr-np.sqrt(float(disc)))/2, float(tr+np.sqrt(float(disc)))/2)

def run(gap, coupling=F(0), eta=F(1,5), eps=F(1,10)):
    # P has a tunable small gap; R is exactly -eta P-eps I plus a PSD slack.
    P=[[F(0),F(0)],[F(0),gap]]
    slack=[[eps,F(0)],[F(0),eps]]
    R=[[-eta*P[i][j]-eps*(i==j)+slack[i][j] for j in range(2)] for i in range(2)]
    # Add an optional rational off-diagonal perturbation while preserving the
    # PSD slack by using slack = vv^T (represented exactly).
    if coupling:
        v=[coupling,coupling]
        slack=[[v[i]*v[j] for j in range(2)] for i in range(2)]
        R=[[-eta*P[i][j]-eps*(i==j)+slack[i][j] for j in range(2)] for i in range(2)]
    H=[[P[i][j]+R[i][j] for j in range(2)] for i in range(2)]
    S=[[R[i][j]+eta*P[i][j]+eps*(i==j) for j in range(2)] for i in range(2)]
    exact_S_min=eig2(S)[0]
    exact_H=eig2(H)[0]
    normR=max(abs(x) for x in eig2(R))
    crude=-normR # Eref=0, P min=0
    certified=-float(eps)
    return dict(gap=str(gap), coupling=str(coupling), eta=str(eta), eps=str(eps),
                S_min=exact_S_min, H_min=exact_H, relative_lower=certified,
                crude_norm_lower=crude, relative_width=exact_H-certified,
                crude_width=exact_H-crude, exact_rational_slack=(S==slack))

def main():
    rows=[]
    for n in [1,2,4,8,16]:
        rows.append(dict(size=n, **run(F(n,10), F(1,10))))
    # Gap closing: relative proof remains valid as gap -> 0; exact H may mix.
    rows += [dict(size=1, **run(F(1,10**k), F(1,20))) for k in [1,2,4,8]]
    # Interacting chain attempt: P is a nearest-neighbour domain-wall parent,
    # while vv^T couples every computational basis state (not disjoint blocks).
    for L in [4,6,8]:
        d=2**L; eta=F(1,5); eps=F(1,10)
        P=np.zeros((d,d),dtype=object)
        for z in range(d):
            bits=[(z>>i)&1 for i in range(L)]
            P[z,z]=sum(bits[i]^bits[(i+1)%L] for i in range(L))
        v=[F(1,100) if z else F(1,20) for z in range(d)]
        H=np.array([[P[i,j]-eta*P[i,j]-eps*(i==j)+v[i]*v[j] for j in range(d)] for i in range(d)],dtype=object)
        # Exact rational lower and rational trial-state upper (u=e_0).
        upper=H[0,0]; exact=float(np.linalg.eigvalsh(np.array(H,dtype=float))[0])
        rows.append(dict(size=L, family='interacting_cycle', exact_ground_diagnostic=exact,
                         relative_lower=-float(eps), rational_upper=float(upper),
                         certified_width=float(upper+eps), premise='S=vvT exact rational',
                         disjoint_blocks=False))
    out='/Users/aidenlippert/Documents/Spectra/results/all_angles_20260913/relative_continuation'
    os.makedirs(out,exist_ok=True)
    with open(out+'/receipt.json','w') as f: json.dump({'units':'abstract','rows':rows},f,indent=2)
    print(json.dumps(rows,indent=2))
if __name__=='__main__': main()
