# Identities and falsification targets

## A. Exact response with collective closure

The user's block congruence is exact. With D>=delta I and delta>0,

    [[K,E*],[E,D]]
      = diag(K-E*E/delta, D-delta I)
        + delta [E/delta, I]* [E/delta, I].

The remaining difficulty is certifying the two diagonal terms without
solving the original problem. This is consistent with established
[Feshbach-Schur analysis](https://arxiv.org/abs/2105.02058); the congruence
itself is not a new theorem.

Consider M bath spins and one central spin, with n=spin-up occupation,
J_+=sum_i s_i^+, J_-=J_+*, and fixed q=n_c+N_b. Let

    H = epsilon n_c + omega N_b + chi n_c N_b
        + g (s_c^+ J_- + s_c^- J_+).

We restrict 1<=q<=M and g!=0. In the central-down/up decomposition,
A0=omega q and D0=epsilon+(omega+chi)(q-1) are scalars, and B=g J_+.
For a trial lower b with d=D0-b>0, choose

    X = (g/d) J_-,    E=0,
    K = (A0-b) I - (g^2/d) J_+J_-.

This is one collective response and one scalar division, with no inverse
on the binomially large eliminated sector. On the bath-q sector, put
c=q(M-q+1), Jmax=M/2. Then

    c I-J_+J_- = Jmax(Jmax+1)I-J^2
               = sum_{i<j} (I-Swap_ij) >= 0.

The last equality follows from Swap_ij=2 s_i.s_j+I/2 and
J^2=3M/4+2 sum_{i<j} s_i.s_j. Each I-Swap is positive because
Swap is a Hermitian involution. Therefore

    K = [(A0-b)-g^2 c/d] I
        + (g^2/d) sum_{i<j}(I-Swap_ij) >=0

whenever (A0-b)d>=g^2 c. This is a universal symbolic sum; the checker
does not expand its M(M-1)/2 factors.

The symmetric Dicke vector d_q=sum_{|S|=q}|S> saturates the collective
bound. The retained/up trial vector (d_q, t J_-d_q) has Rayleigh quotient

    U(t) = [A0 + t^2 c D0 + 2gtc] / [1+t^2 c].

Rational t gives an exact upper. Its combinatorial normalization cancels;
the checker need not construct binomial coefficients or amplitudes. The
lowest root of (A0-b)(D0-b)=g^2 c is the exact ground energy in this sector.
Rational bisection and the displayed trial give a rigorous interval.

This is an explicitly solvable homogeneous central-spin family, related
to established [XXZ central-spin solutions](https://arxiv.org/abs/1810.03012),
not a new molecular algorithm. At fixed precision its scalar checker uses
O(log M)-bit parameter arithmetic, while the expanded Hamiltonian contains
M couplings and the sector contains binomial(M,q)+binomial(M,q-1) states.
The compression statement assumes the Hamiltonian is supplied by the
homogeneous formula. Recognizing an arbitrary expanded input costs at least
its input length. Applying J_- to an arbitrary explicit vector is not cheap.

Falsification targets on H6: (1) Does its pattern algebra close in a tiny
commuting collective algebra? Generate exact commutators of its ten spatial
6x6 patterns. Full rank 35 proves sl(6,Q) closure and refutes that hypothesis.
It does not rule out every noncommutative response. (2) For the concrete split
by occupation of spin orbital zero, is D scalar? Two unequal determinant
diagonals inside D suffice to disprove it. Neither test bounds a general
response residual or supplies an eliminated-sector gap certificate.

Fresh controls: test the exact response at new M and detunings, including
large M without enumeration. At M=6 deliberately make one coupling
inhomogeneous. The scalar-D response still exists, but the symmetric
Casimir bound need not apply to the new weighted B. Never apply the
homogeneous checker to that changed Hamiltonian.

A precise possible extension beyond scalar D is an intertwining identity.
Suppose Q_e L=L(Q_r+kappa I), A=a(Q_r)-bI, D=f(Q_e)-bI, B*=L,
and phi(Q_r)=f(Q_r+kappa)-bI is strictly positive. Functional calculus gives
D L=L phi, so X=L phi^-1 has E=0. Taking the adjoint of the intertwining
identity proves that L*L commutes with Q_r. Consequently

    K = a(Q_r)-bI - L*L phi(Q_r)^-1.

If L*L<=c(Q_r), the remaining certificate is a scalar functional inequality
a(x)-b-c(x)/phi(x)>=0 on the allowed charge spectrum, together with the
independent bound D>=delta I. Successive intertwiners add their shifts.
This derivation does not establish that a short molecular L, a tractable
charge spectrum, or the required norm and gap certificates exist. Multiple
branches can proliferate under composition. No checker for this broader
class was built because the current diagnostic establishes no useful
molecular instance. The central-spin control realizes only its scalar case.

## B. An existing-moment consistency law outside the old cone

For any alternating three-form c on the one-particle space, define

    C(c)=sum_{i<j<k} c_ijk a_i a_j a_k,
    P(c)=C(c)* C(c)+C(c) C(c)* >=0.

For any vector psi, its expectation is ||C psi||^2+||C* psi||^2.
Moving three creation operators past three annihilation operators takes
nine swaps. Thus the uncontracted sixth-degree terms have opposite signs
in the two products. All surviving terms have degrees 0, 2 or 4. No
three-particle density matrix or guessed sixth-degree moments are needed.
This is the known T1 condition; see
[Zhao, Braams, Fukuda, Overton and Percus](https://optimization-online.org/2003/10/760/).

Equivalently, the matrix T1_y[I,J]=y({a_I*,a_J}) must be positive for
every physical degree-four moment functional y. This is a parametric law
for every three-form, every orbital basis, and every fermionic state, not
a theorem specific to the fitted coefficients below. The old charge-one
patterns do not implement the complete T1 class. Strict non-implication is
established by the *full old-family audit plus an exact negative evaluation*,
not merely by counting or renaming generators.

The parent's checker defines a sparse rational functional on all balanced
canonical words of degrees 0, 2, 4, with omitted values zero. We rebuild this
domain explicitly (1+m^2+binomial(m,2)^2 words). A candidate outside it is
refused. We do not extend that convention to unspecified higher moments.

If y(P)<0 then no positive many-body state, and no admissible higher-moment
extension of y, can exist: any such extension would assign nonnegative
expectations to both squares while their sum is already fixed and negative.

Robustness test fixed before further evaluation: with tau the exact uniform
fixed-N trace functional, y_t=(1-t)y+t tau remains feasible for every old
rule by convexity. Evaluate the frozen C at t=1/1000,1/200,1/100,1/50.
This gives fresh controlled points in a structural class, not independent
molecular transfer. Also test the identity on fresh rational physical states
with an independently implemented fermionic action. Charge their enumeration.
No improved energy bound or scalable search is inferred from a separator.
