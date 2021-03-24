# SIR model with testing and isolation mechanisms, (by Ali Gharouni)
# non zero p_S was assumed from the begining in the equations

restart:
with(LinearAlgebra):
with(linalg):
with(plots):
#sc:=tau*rho*N_0/(tau*(W_s*S_u+W_i*I_u+W_r*R_u)+rho*N_0);# scale
 # note that in the Rmd file, t is replaced by rho

sc:=rho*N_0/(W_s*S_u+W_i*I_u+W_r*R_u);
F_S:=sc*W_s; F_I:=sc*W_i; F_R:=sc*W_r;
# the FOI

Lambda:= beta*(I_u+eta[w]*(I_n+I_p)+eta[c]*I_c)/N_0;
# The equations of the model:

eq1:= -Lambda*S_u - F_S*S_u+ omega*S_n;
eq2:= -Lambda*S_n + (1-p_S)*F_S*S_u - omega*S_n;
eq21:= -Lambda*S_p + p_S*F_S*S_u - omega*S_p;
eq22:= -Lambda*S_c + omega*S_p;
eq3:= Lambda*S_u + omega*I_n - F_I*I_u - gamma*I_u; 
eq4:= Lambda*S_n + (1-p_I)*F_I*I_u - omega*I_n -gamma*I_n;
eq5:= Lambda*S_p + p_I*F_I*I_u - omega*I_p -gamma*I_p;
eq6:= Lambda*S_c + omega*I_p - gamma*I_c;
eq7:= gamma*I_u + omega*R_n - F_R*R_u;
eq8:= gamma*I_n + (1-p_R)*F_R*R_u - omega*R_n;
eq9:= p_R*F_R*R_u + gamma*I_p - omega*R_p;
eq10:= gamma*I_c + omega*R_p;
# Sanity Check:
simplify(eq1+eq2+eq21+eq22+eq3+eq4+eq5+eq6+eq7+eq8+eq9+eq10);


# Stability of Disease Free Equilibrium When p_S=0 is assumed (that is the prob that a S individual  be tested positive, thus eq21 and eq22 are not included here)
with(linalg):
with(LinearAlgebra):
with(VectorCalculus): # We need this library to calculate the Jacobian matrix.

with(ArrayTools):
# The bifurcation parameter is assumed to be the testing intensity, rho
G[1]:= (S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau) -> eq1:
G[2]:= (S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau) -> eq2:
G[3]:= (S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau) -> eq3:
G[4]:= (S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau) -> eq4:
G[5]:= (S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau) -> eq5:
G[6]:= (S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau) -> eq6:
G[7]:= (S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau) -> eq7:
G[8]:= (S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau) -> eq8:
G[9]:= (S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau) -> eq9:
G[10]:=(S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau) -> eq10:
J:= Jacobian([ 
G[1](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[2](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[3](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[4](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[5](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[6](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[7](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[8](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[9](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[10](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau)
],
[S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c]):
J;
# Difine the Jacobian as a function of all params:

JJ := unapply(J,S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau):
                

# Jacobian at DFE, ie when all I and R compartments are = 0 ;

J_eval:=JJ(Su,Sn,0,0,0,0,0,0,0,0,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau);
# algsubs(S_n+S_u=N_0,JJ(S_u,S_n,0,0,0,0,0,0,0,0,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c]));
# Eigenvalues(J_eval);

# What is the eigenvalues of the last 4 equations, the recovery?
J_subR:=Jacobian([ 
G[7](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[8](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau)
],[R_u,R_n]);
JJ_subR := unapply(J_subR,S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau):
JJ_subR(S_u,S_n,0,0,0,0,0,0,0,0,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau);
# Both eigenvalues are negative.

Eigenvalues(JJ_subR(S_u,S_n,0,0,0,0,0,0,0,0,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau));
J_subI:= Jacobian([ 
G[3](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[4](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[5](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[6](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau)
],[I_u,I_n,I_p,I_c]);


JJ_subI := unapply(J_subI,S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau):
simplify(JJ_subI(S_u,S_n,0,0,0,0,0,0,0,0,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau));
#Eigenvalues(JJ_subI(Su_dfe,Sn_dfe,0,0,0,0,0,0,0,0,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c]));
#CharacteristicPolynomial(JJ_subI(Su,Sn,0,0,0,0,0,0,0,0,rho,omega,beta,gamma,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),lambda);
# Let's define Matrix JJ_subI as simple as possible and see how the eigenvalues are:
JJ_subI2:=Matrix([[a1,a2,a3,a4],[b1,b2,b3,b4],[p_I*F_I,0,-omega-gamma,0],[0,0,omega,-gamma]]);
CharacteristicPolynomial(JJ_subI2,lambda);
# the I_u and I_n gives another decoppled part of the whole Jacobian Matrix
J_sub2:= Jacobian([ 
G[3](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_S,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau),
G[4](S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_S,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau)
],[I_u,I_n]);

#JJ_sub2 := unapply(J_sub2,S_u,S_n,I_u,I_n,I_p,I_c,R_u,R_n,R_p,R_c,rho,omega,beta,gamma,p_S,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c],tau);
#JJ_sub2(S_u,S_n,0,0,0,0,0,0,0,0,rho,omega,beta,gamma,p_S,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c]);
#(collect(CharacteristicPolynomial(JJ_sub2(S_u,S_n,0,0,0,0,0,0,0,0,rho,omega,beta,gamma,p_S,p_I,p_R,W_s,W_i,W_r,N_0,eta[w],eta[c]),lambda),lambda));


# SIR with isolation parameters eta[w] and eta[c]
restart:
with(LinearAlgebra):
with(linalg):
with(plots):
with(DEtools):
F:= Matrix([
[beta*S_u/N_0, beta*eta[w]*S_u/N_0, beta*eta[w]*S_u/N_0, beta*eta[c]*S_u/N_0],
[beta*S_n/N_0, beta*eta[w]*S_n/N_0, beta*eta[w]*S_n/N_0, beta*eta[c]*S_n/N_0],
[beta*S_p/N_0, beta*eta[w]*S_p/N_0, beta*eta[w]*S_p/N_0, beta*eta[c]*S_p/N_0],
[beta*S_c/N_0, beta*eta[w]*S_c/N_0, beta*eta[w]*S_c/N_0, beta*eta[c]*S_c/N_0]
]);
V:= Matrix([
[F_I+gamma, -omega, 0,0],
[-(1-p_I)*F_I, omega+gamma,0,0],
[-p_I*F_I, 0,omega+gamma,0],
[0,0,-omega,gamma]
]);
# V^-1

Vi:=MatrixInverse(V);
# Assumption: p_S=0, thus S_p=S_c=0

G:=subs([S_p=0,S_c=0],simplify(Multiply(F,Vi)));
# Charactristic poly in x:

collect(simplify(CharacteristicPolynomial(G,x)),x);
# Method1 (symbolically): R0 is given by

sol:=solve(CharacteristicPolynomial(G,x),x); 
# Method 2 (by hand) check if the formula I derived for R0 is valid
A:=gamma*(omega+gamma+eta[w]*F_I)+omega*eta[c]*p_I*F_I;
B:=(gamma*(omega+gamma)*(omega+eta[w]*(gamma+F_I))+omega*p_I*F_I*(gamma*eta[w]+omega*eta[c]))/(omega+gamma);
C:=beta/((gamma*(omega+gamma)+F_I*(gamma+omega*p_I))*gamma*N_0);

collect(B,F_I);
R0:=simplify((S_u*A+S_n*B)*C);
# Cheching if the symbolic computation of R0 and by hand one are the same
simplify(sol[4]-R0);
# Derivative of R0 wrt parameters
with(plots):
diff(collect(R0,omega), omega);
# Define some key parameters:

N0:=1000;
sc := rho*N_0/(W_s*S_u+W_i*I_u+W_r*R_u);
#At Disease Free Equilibrium:
Sn_dfe:=subs(p_S=0,(rho*(1-p_S)*N_0)/omega);
Su_dfe:=N_0-Sn_dfe;



Fs_dfe := simplify(rho*N_0/Su_dfe);
Fi_dfe := simplify(rho*N_0*W_i/(W_s*Su_dfe));
params_ini := [ beta=1,gamma=1/3, omega=0.25, rho=0.01, 
             W_s=0.01, W_i=1, W_r=0.01,
             p_S=0, p_I=0.5, p_R=0.5, eta[w]=0.02, eta[c]=0.019, 
             N_0=N0];
params_dfe := [ params_ini[], 
             S_u=eval(Su_dfe, params_ini),S_n=eval(Sn_dfe, params_ini),I_u=1,R_u=0,
             F_S=eval(Fs_dfe, params_ini), F_I=eval(Fi_dfe, params_ini)];
# Factoring out eta 1 and 2, Consider a parametrization s=eta[w]/eta[c] note that s>=1 is our interest
algsubs(S_n+S_u=N_0,collect(simplify(R0),[eta[w],eta[c]]));
Sn_dfe;
Su_dfe;
Fi_dfe;
subs([S_u=Su_dfe,S_n=Sn_dfe,F_I=Fi_dfe],R0) assuming W_s>0, tau>0, omega>0;
# Tylor expansion of R0 at t=0
tl_t:=simplify(taylor(subs([S_u=Su_dfe,S_n=Sn_dfe,F_I=Fi_dfe],R0),rho=0,3)) assuming W_s>0, tau>0, omega>0; # Notice that the degree 2 for taylor series doesn't work !
# derivative of taylor expansion of R0 wrt omega
d_omega:=diff(simplify(subs([S_u=Su_dfe,S_n=Sn_dfe,F_I=Fi_dfe],R0)),omega); # Mess!
simplify(diff(convert(tl_t,polynom),omega));
collect(simplify(-1/(beta*rho)*collect(numer(simplify(diff(convert(tl_t,polynom),omega))),omega)),omega);

#implicitplot3d(subs(params,R0)=1, omega=0.1..1,eta[c]=0.1..1,rho=0..0.01,axes=boxed,orientation=[-92,71,2]);

# I want to show that diff (BCSn*) is negative wrt rho
collect(simplify(subs(F_I=Fi_dfe,B*C*Sn_dfe)),rho);
collect(simplify(diff(simplify(subs(F_I=Fi_dfe,B*C*Sn_dfe)),rho)),rho);
# collect t and F_I

BCSn:=collect(simplify(subs(B*C*Sn_dfe)),[rho,F_I]);
# let W=W_I/W_S, note W>=1

Fi_dfe2:=algsubs(W_i/W_s=W,Fi_dfe);
# subs in the Fi_dfe into BCS_n and try to simplify the terms and then diff wrt rho;
# simplify the Numerator of BCSn 
(collect(simplify(subs(F_I=Fi_dfe2,numer(BCSn))),[rho]));
n1:=collect(numer(%),rho);
# understanding the coefficients of n1, decompositions;
factor(coeff(n1,rho));
factor(coeff(n1,rho^2));
collect(coeff(n1,rho^2),-1);

# simplify the Denominator of BCSn 

collect(simplify(subs(F_I=Fi_dfe2,denom(BCSn))),[rho]);
d1:=collect(numer(%),rho);
# root of d1:
d1_root:=solve(d1,rho);
# check if d1_root can be factored in the numerator?
simplify(subs(rho=d1_root,n1));
# conclude that d1_root canno be factored in the numerator. Note the evaluated expression is positive.
d1_simp:=collect(simplify(d1/(gamma*omega*(omega+gamma))),rho);



# DIFF PART:
# seperatly diff wrt t and apply qoutien rule just in the numerator;
diff1:=collect(simplify(diff(n1,rho)*d1_simp-diff(d1_simp,rho)*n1),rho);
collect(subs([p_I=1],factor(diff1)),rho);






# 















