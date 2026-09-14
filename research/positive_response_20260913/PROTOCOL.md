# Nonzero response error and a molecular gap obstruction

Preserve the preceding 718-file manifest. No molecular SDP solve, adaptive
direction enrichment, or new fixture/FCI generation is planned in this pass.
The repository already has frustration-free controls and enumerated Schur
responses. Reuse their mathematics and label the present control's scope.

First derive the following exact alternative to paying ||E||^2/delta.
If a reference H0-b admits response X0 with E0=0, K0>=0 and D0>0, and
H=H0+diag(V0,V1), V0,V1>=0, then

 T*(H-b)T = diag(K0,D0)+diag(V0,0)+[-X0,I]* V1 [-X0,I].

This certifies the actual nonzero E=-V1 X0 collectively. It does not require
V1 to commute with X0, nor small ||E||. Test it on a central spin plus a
conditional exchange bath V_s=lambda_s sum_i n_(i+2)(I-Swap_(i,i+1)), with
periodic bath indices. Each local term is positive and kills every symmetric
Dicke state. Thus the preceding scalar reference supplies a matched upper
too. The model is deliberately frustration-free in its added interactions;
it is a proof-rule control, not a generic molecular mechanism. Explicitly
test noncommutation, nonzero residual and failure of the scalar error test.

Second derive the molecular limitation quantitatively. For N=s electrons
in s spatial orbitals, the existing orbital/spin Casimir identity becomes

 C_orb=sum_pq (E_pq-delta_pq N/s)* (E_pq-delta_pq N/s)
       =2[Smax(Smax+1)-S^2], Smax=s/2.

If the retained density patterns generate sl(s), their common kernel is
exactly the maximum-spin multiplet. Compare its energy with an already
certified physical upper. This tests whether treating those interactions
as a bulk with a common zero can describe the molecule's low energy.

Try one explicit quantitative gap certificate. A bracket-generated basis
gives ||Q(A_j)psi||^2<=k_j <R>, R=sum lambda_i Q(L_i)^2/2. Use elementary
operator norm bounds and exact rational linear reconstruction to bound the
full Casimir by K R. Then R>=C_orb/K, and R>=2s/K on any subspace with a
doubly occupied orbital (which excludes Smax). Record all norm constants,
coefficient sizes and construction/replay costs. Test the resulting cheap
gap for each of the six double-occupancy complements at the frozen H6
target. Failure of this sufficient bound is not proof that D lacks a gap.
Do not optimize a new molecular certificate family if the bound is weak.

As a consistency check, project the frozen T1 three-removal operator onto
its spin-one-half component by the exact SU(2) operator Casimir. This is a
symmetry projection, without coefficient fitting. Verify spin covariance,
sixth-degree cancellation and the negative moment using existing moments
only. Determine its action on the maximum-spin common kernel. The T1 law
and Casimir identities are established mathematics; report new instance
results and limitations separately from those identities.

Fresh controls use new positive-bulk strengths and noncommuting conditional
exchange instances. The H8 stored fixture is held out from the H6 algebra
calculation for a kernel/energy diagnostic. Charge all small basis
enumerations. Accepting paths use exact rational arithmetic and Python -S.
