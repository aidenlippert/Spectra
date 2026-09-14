# A positive block can absorb a large response error

Let H0-b=[[A_ref,B],[B*,D_ref]], let X0=D_ref^-1 B* denote an already available
reference response (the implementation uses a scalar D0), and put
K0=A_ref-BX0>=0. The inverse notation here specifies the identity; it is not
an instruction to evaluate a many-body inverse. With
T=[[I,0],[-X0,I]], and any V0,V1>=0,

 T*(H0+diag(V0,V1)-b)T
 = diag(K0,D_ref)+diag(V0,0)+[-X0,I]* V1 [-X0,I] >=0.

Direct expansion gives E=-V1X0 and K=K0+V0+X0*V1X0. A global error test
would instead require K>=X0*V1^2X0/delta. That sufficient test can fail as
V1 grows, even though the displayed block is positive for every V1>=0.
The benefit is preserving the exact relationship among the error terms.
This is ordinary congruence/SOS positivity, not a newly claimed theorem.

For the homogeneous central-spin reference from the previous pass,
H0=epsilon*n_c+omega*N_b+chi*n_c*N_b+g*(sigma_c^+ J_-+sigma_c^- J_+),
fix n_c+N_b=q. Its unshifted block energies are A_energy=omega*q and
D_energy=epsilon+(omega+chi)*(q-1). The response is X0=(g/d)J_- with
d=D_energy-b>0. Define the bath operator

 L=sum_i n_(i+2) (I-Swap_(i,i+1)),  indices modulo M.

The three sites are distinct for M>=3. The occupation projector commutes
with the swap; hence each summand is positive. Equivalently it is
F_i*F_i/2 with F_i=n_(i+2)(I-Swap_(i,i+1)). Every symmetric Dicke vector
is killed by I-Swap, so L kills the retained and eliminated reference
trial vectors. Set V0=lambda0 L, V1=lambda1 L, with nonnegative rational
strengths. D=dI+V1 has the explicit independent lower dI. It is generally
nonscalar, [L,J_-] is generally nonzero, and E is generally nonzero. The
reference lower and rational trial upper still enclose the actual ground
energy. Checking a formula-specified instance needs no many-body expansion;
reading an arbitrary list of local terms would still cost its input length.

The shared-kernel property is a strong restriction. Existing repository
controls already exploit it; see research/certificate_scaling/
low_rank_structural_routes.md. The additional control here tests a
nonzero elimination residual rather than claiming large-system novelty.

# The retained molecular bulk cannot vanish at low energy

Write E_pq=sum_sigma a_p,sigma* a_q,sigma and R=sum_i lambda_i Q(L_i)^2/2.
All lambda_i are strictly positive and every L_i is traceless Hermitian.
If R psi=0, all Q(L_i)psi=0. Their commutators also annihilate psi because
[Q(A),Q(B)]=Q([A,B]). If they generate sl(s), every traceless Q(A) kills
psi. The existing exact CAR identity is

 C_orb=sum_pq (E_pq-delta_pq N/s)* (E_pq-delta_pq N/s)
      =(s+2)N-(1/2+1/s)N^2-2S^2.

On N=s this equals 2[Smax(Smax+1)-S^2], Smax=s/2. Its kernel is the
maximum-spin multiplet. It has exactly one electron per spatial orbital;
its orbital wavefunction is the one-dimensional determinant representation.
Conversely that multiplet is annihilated by every traceless spatial Q.
Thus the kernel of R has dimension s+1 and cannot contain a doubly
occupied orbital. On this kernel the spin-independent retained Hamiltonian
has the scalar energy c_one+trace(t_spatial).

An exact energy above a known physical upper excludes this kernel from the
ground space. It does not prove that every shifted/frustrated SOS
representation fails, or that a different positive bulk has no useful kernel.

# A quantitative commutator gap, with explicit costs

For any spatial matrix A and N particles,
||Q(A)||<=N||A||_2<=N max(||A||_1,||A||_infinity)=M(A).
For normalized generator A_i=L_i/r_i, r_i the last matrix norm bound,

 ||Q(A_i)psi||^2 <= [2/(lambda_i r_i^2)] <psi,R psi>.

If A=[A_i,A_j]/r, then

 k(A)=2[M(A_i)^2 k(A_j)+M(A_j)^2 k(A_i)]/r^2

is valid by the triangle and squared-norm inequalities. No state is
enumerated in this recursion. Choose an independent bracket basis, and
reconstruct the canonical traceless matrices B_t in it: B_t=sum_j c_tj A_j.
For rational u_j>=sqrt(k_j),

 sum_t ||Q(B_t)psi||^2
 <= sum_t (sum_j |c_tj|u_j)^2 <R> = K <R>.

Choose the off-diagonal matrix units and diagonal E_ii-E_(s-1,s-1).
Their summed operator squares dominate C_orb because the inverse diagonal
Gram is I-11^T/s<=I. Therefore R>=C_orb/K. Every state with a specified
doubly occupied orbital has S<=s/2-1, so C_orb>=2s there and R>=2s/K.
This is a nonenumerating positive gap for that interaction bulk. Whether
it is quantitatively useful for the *full* D remains a numerical question:
the one-body term and the fixed tail must be included before subtracting b.

The checker rebuilds the bracket identities, norm constants and linear
reconstructions using rationals. Triangle bounds and ill-conditioned
reconstruction may make K large. Its cost and weakness must be reported.

For an actual molecular eliminated space choose Q_p=n_(p,up)n_(p,down).
The compressed one-body operator on this space is
c_one+2t_pp+dGamma(t_without_p), with s-2 remaining electrons.
Let rho_p be the maximum absolute off-diagonal row sum of t_without_p.
Its off-diagonal part is bounded below by -rho_p I, hence

 Q_p H Q_p >= [c_one+2t_pp
              +2 sum_(lowest (s-2)/2 diagonal entries) t_ii
              -(s-2)rho_p+ell+2s/K] Q_p.

Here ell is the already certified lower tail shift. This elementary
bound can also be used with 2s/K replaced by zero. The target is fixed at
b=U-1/625 Ha, using the frozen physical upper U. A positive difference
between this lower bound and b proves the required D gap; it does not
bound the retained Schur remainder or construct a response. Negative
differences are inconclusive about the actual gap. The excluded sector
has binomial(2s-2,s-2) states; the formula is counted without enumerating it.

# Spin projection of the existing consistency defect

For the frozen C, [S_z,C]=-C/2. The adjoint SU(2) Casimir

 J_ad^2(C)=[S_z,[S_z,C]]
           +([S_+,[S_-,C]]+[S_-,[S_+,C]])/2

has eigenvalues 3/4 and 15/4 on three spin-one-half annihilators. Hence
C_half=(15 C/4-J_ad^2(C))/3 is the exact spin-one-half component.
It obeys [S_-,C_half]=0; the three-particle/hole anticommutator remains
quartic and positive. Its negative parent moment must be checked anew.

C_half and its adjoint annihilate the N=s, S=s/2 multiplet: removing or
adding three electrons leaves maximum spin (s-3)/2, which cannot be
reached from s/2 by a rank-one-half operator. An independent finite action
control also tests this statement. The projected T1 defect therefore
concerns the lower-spin correlations that the common-zero bulk omits.
It is still an instance of the established T1 condition, not a new law.

The spin partner C_plus=[S_+,C_half] completes the doublet. Therefore

 P_scalar={C_half*,C_half}+{C_plus*,C_plus} >=0

commutes with S_z and S_+/-; exact CAR expansion checks these identities.
This removes spin-axis dependence from the separator. In particular its
negative expectation is unchanged by any global spin rotation of the
functional. Old-family feasibility of those rotated functionals is not
assumed here. Convex mixtures with the uniform fixed-N trace do remain
old-family feasible, using the preceding exact parent audit. For
y_t=(1-t)y+t*tau, the strict rejection threshold is
-y(P_scalar)/(tau(P_scalar)-y(P_scalar)), computed exactly.
