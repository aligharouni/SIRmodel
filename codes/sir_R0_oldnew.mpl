# SIR model with testing and isolation mechanisms, (by Ali Gharouni)
# Note on 27 Apr 2021, we decided to consider a reduction scaler in S_n to I_n

# 
# SIR old notation; R0 with isolation parameters eta[w] and eta[c]
restart:
with(LinearAlgebra):
with(linalg):
# Forming the next generation matrix G=FV^-1
F:= Matrix([
[beta*S_u/N_0, beta*eta[w]*S_u/N_0, beta*eta[w]*S_u/N_0, beta*eta[c]*S_u/N_0],
[beta*S_n*eta[w]/N_0, beta*eta[w]^2*S_n/N_0, beta*eta[w]^2*S_n/N_0, beta*eta[c]*eta[w]*S_n/N_0],
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
# Method 2 (component equation; derived by hand) check if the formula I derived for R0 is valid
A:=gamma*(omega+gamma+eta[w]*F_I)+omega*eta[c]*p_I*F_I;
B:=(gamma*(omega+gamma)*(omega+eta[w]*(gamma+F_I))+omega*p_I*F_I*(gamma*eta[w]+omega*eta[c]))/(omega+gamma);
C:=(N_0*(gamma*(omega+gamma)+F_I*(gamma+omega*p_I)));
# R0; Basic reproduction number, note S_n is scaled by eta_w now
R0:=simplify(beta/gamma*(S_u*A+S_n*eta[w]*B)/C);
# Cheching if the symbolic computation of R0 and by hand one are the same
simplify(sol[4]-R0);
subs(N_0=1,simplify(sol[4]-R0));
# simplify and factor etas' from the numerator of R0
collect((S_u*A+S_n*eta[w]*B),[eta[w]]);


# 

# 
# Date 27 Apr 2021, isolation parameters 
# The symbols are based on Appendix of manuscript.
# this is C in the appendix

C0:=((omega+gamma)*(gamma*omega+gamma^2+F_I*gamma+omega*p_I*F_I));
C1:=(omega+gamma)*(theta[w]*gamma+theta[c]*omega*p_I)*F_I;
C2:=(omega+gamma)*gamma^2*theta[w]+(theta[w]-theta[c])*gamma*omega*p_I*F_I;
collect(simplify(C1+C2),[gamma]);
C2_:=collect(simplify(C1+C2),[theta[w],theta[c]]);



Delta4:=C1/C0*S_u+(C2_/C0*(1-theta[w])+theta[w])*S_n;
R04:=beta/gamma*(1-Delta4);
newR0:=simplify( algsubs(S_u+S_n=1,R04) );
R0_:=simplify(algsubs(eta[c]=1-theta[c],algsubs(eta[w]=1-theta[w],R0)));
oldR0:=subs(N_0=1,algsubs(S_u+S_n=N_0,R0_));
simplify(oldR0-newR0);
# IN the new Delta4
simplify(Delta4);
C0;
# note that the negative will b etaken care of when using numer function on Delta4
numer(simplify(Delta4));
# For some strange reason S_u+S_n=1 was not applied! So we reinforce it here:
algsubs(S_u+S_n=1,numer(simplify(Delta4)));
# Verifying numerator of Delta4 wrt theta:
Delta4_theta:=collect(algsubs(S_u+S_n=1,numer(simplify(Delta4))),[theta[w],theta[c]]);
coef2:=collect(factor(coeff(Delta4_theta,theta[w],2)),[F_I]);
coef1:=collect(coeff(Delta4_theta,theta[w],1),[theta[c],S_u,F_I]);
# Delta4 is face down quadratic in theta_w, with one root 0 and 1 root positive, I want to see if the vertex of this parabola lies beyond theta_w=1, if so there is a range for theta_w st diff Delta4 wrt theta_w becomes negative so more theta_w means more R0!!? PUZZELED! 
# one trick is to consider diff of numer Delta4 wrt theta_w and evaluate it at theta_w=1, to see if it is pos or neg? it gives us a sense where theta_w=1 is located wrt vertex of parabola in theta_w...
collect(simplify(2*1*coef2+coef1),[S_u]);










# Date 5 May 2021, chnages in R0 wrt testing intensity rho
# Note that w=w_I/w_S, here for brevity I use w.
F_Ihat:=rho*omega/(omega-rho)*w;
F_Shat:=rho*omega/(omega-rho); # this is phi
Delta4;
# replace F_I and S_n and S_u at DFE in Delta:  
Delta4_simp:=algsubs(F_I=F_Ihat,algsubs(S_n=rho/omega,algsubs(S_u=1-S_n,Delta4) ));
# Reparametrization of rho bt phi:
Delta4_repar:=simplify(subs(rho=omega*phi/(omega+phi),Delta4_simp));
collect(numer(Delta4_repar),phi);
# Simplifying K1:
K1:=factor(coeff(collect(numer(Delta4_repar),phi),phi,2));
K1/w;
# which factorization is better in K1/w?
collect(K1/w,[theta[c],theta[w]]);
collect(K1/w,[omega,gamma]);

K2:=factor(coeff(collect(numer(Delta4_repar),phi),phi,1));
K3:=factor(denom(Delta4_repar));
#simplify(diff(Delta4_repar,phi));
dDelta4phi:=collect(simplify(diff(Delta4_repar,phi)),phi);
# Just look at the numer with omega factored
dDelta4phi_numer:=collect(numer(dDelta4phi)/omega,phi);
a3:=factor(coeff(dDelta4phi_numer,phi,2));
# Can we factor a3 more?
collect(a3/w,[gamma,omega,p_I,theta[c],theta[w]]);
#collect(a3/(w),[w,omega,gamma]);
# check a3<0 if theta_w=0,theta_c=1 (ie the top right corner of panel b)
simplify(subs(theta[w]=0,theta[c]=1,a3/w));

#latex(%);
# Checking in the case of random testing, i.e., w_I=w_S which is equivalent ti w=1, results in a3>0? Yes, just look at the coef of gamma omega^2 and vary theta_c and p_I from 0 to 1. Realize that in the 4 extreme cases the coef is positive (see page 33 of my notes)
collect(subs(w=1,a3),[omega,gamma]);
# check a3 when w=1 and p_I=1:
collect(subs([w=1,p_I=1],a3),[omega,gamma]);
b3:=factor(coeff(dDelta4phi_numer,phi,1));
# Note that b3 is positive
collect(b3/(2*gamma*(omega+gamma)*w),[theta[c],theta[w]]);
latex(%);
# Note that c3 is positive
c3:=factor(coeff(dDelta4phi_numer,phi,0));
# what happens in theta_w=0 theta_c=1?
factor(subs(theta[c]=1,subs(theta[w]=0,a3)));
# Linearization of  delta at rho=0, where R0=beta/gamma(1-Delta)
Delta4_simp;
tl_t:=simplify(taylor(Delta4_simp,rho=0,2)) assuming w>0, rho>0, omega>0;
collect(numer(convert(tl_t,polynom)),[rho,omega,gamma]);
collect(simplify(diff(convert(tl_t,polynom),omega)),omega);







# R0 wrt w=w_I/w_S
Delta4;
Delta4_simp;
dD4dw:=simplify(diff(Delta4_simp,w));
factor(numer(dD4dw)/(gamma*rho));
# Main point here is: The latter expression is always positive so D Delta/D w>0 so the more targeted testing the lower R0

# What is the denominator?
denom(dD4dw);
factor(-omega^2*gamma+omega*gamma*rho-gamma*rho*omega*w-omega*gamma^2+gamma^2*rho-omega^2*p_I*rho*w);
latex(%);






















