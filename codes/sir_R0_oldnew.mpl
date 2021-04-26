# SIR model with testing and isolation mechanisms, (by Ali Gharouni)

# SIR old notation; R0 with isolation parameters eta[w] and eta[c]
restart:
with(LinearAlgebra):
with(linalg):
# Forming the next generation matrix G=FV^-1
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
# Method 2 (component equation; derived by hand) check if the formula I derived for R0 is valid
A:=gamma*(omega+gamma+eta[w]*F_I)+omega*eta[c]*p_I*F_I;
B:=(gamma*(omega+gamma)*(omega+eta[w]*(gamma+F_I))+omega*p_I*F_I*(gamma*eta[w]+omega*eta[c]))/(omega+gamma);
C:=1/(N_0*(gamma*(omega+gamma)+F_I*(gamma+omega*p_I)));
# R0; Basic reproduction number
R0:=simplify(beta/gamma*(S_u*A+S_n*B)*C);
# Cheching if the symbolic computation of R0 and by hand one are the same
simplify(sol[4]-R0);
subs(N_0=1,simplify(sol[4]-R0));



# SIR new notation; R0 = beta/gamma(1-Delta) with theta symbols
# Components of R0

A_:=collect(S_n*omega*gamma^2+gamma^2*F_I+S_n*omega*p_I*F_I*gamma+S_n*gamma^3+omega*gamma*F_I,[S_n]);
B_:=omega*p_I*F_I* (gamma*S_u + omega );
C_:=((omega+gamma)*(gamma*omega+gamma^2+F_I*gamma+omega*p_I*F_I));
A_ * theta[w] + B_ * theta[c];
Delta:=%/C_;
R0_:= simplify( algsubs(S_u+S_n=1,beta/gamma * (1-Delta)) );
simplify(algsubs(theta[c]=1-eta[c],algsubs(theta[w]=1-eta[w],R0_)));
# Check between old and new notation; checked!
simplify(%-subs(N_0=1,algsubs(S_u+S_n=N_0,R0)));


# 
C_;
Visub:=Vi[..,1..2];
Visub2:=Matrix([
[gamma*(omega+gamma)^2,gamma*omega*(omega+gamma)],
[gamma*(omega+gamma)*(1-p_I)*F_I,gamma*(omega+gamma)*(F_I+gamma)],
[gamma*(omega+gamma)*p_I*F_I,gamma*omega*p_I*F_I],
[omega*(omega+gamma)*p_I*F_I,omega^2*p_I*F_I]
]);
1/(gamma*C_)*Visub2;
simplify(Visub- 1/(gamma*C_)*Visub2);
F1:=Matrix([[S_u],[S_n],[0],[0]]);
F2:=Matrix([1,1-theta[w],1-theta[w],1-theta[c]]);
F0:= Multiply(F1,F2);
beta*F0;

simplify(Multiply(F2,Visub2));






# New representation of R0 when Delta is written in terms of rho and omega:
# R0=beta/gamma(1-Delta); 
# old representation with factors of theta.
Delta;
C_1:= (omega+gamma)*(theta[w]*gamma+theta[c]*omega*p_I)*F_I;
C_2:=gamma^2*theta[w]*(omega+gamma)+gamma*omega*(theta[w]-theta[c])*p_I*F_I; 
# is Delta=(K_1+K_2)S_n/C?
C_; 
N_0:=1;
F_Ihat:=rho*omega/(omega-rho)*w;
F_Shat:=rho*omega/(omega-rho);
# Looking at matrix G_11 and breaking down the 2 columns of FV^-1, ie C-C1 and C-C1-C2 and see how they change wrt rho?
simplify(algsubs(F_I=F_Ihat,(C_-C_1)));


collect(algsubs(F_I=F_Ihat,(C_-C_1)),rho);





# A new representation of Delta to simplify the diff Delta wrt rho and omega I hope.
Delta2:=(C_1+C_2*S_n)/C_;
simplify(Delta2);
simplify(algsubs(S_n+S_u=N_0,Delta));
# We only need to check the numerators of Delta and Delta2
numer(Delta2);
collect(numer(Delta2),S_n);
numer(Delta);
collect(simplify(algsubs(S_u=1-S_n,numer(Delta))), S_n);
# check if the two representations match?
simplify(algsubs(S_n+S_u=1,numer(Delta)-numer(Delta2)));
simplify(algsubs(S_n+S_u=1,Delta-Delta2));
Delta2;
F_Ihat;
# let w:=w_I/w_S, and substitude F_I (this is \hat F_I in the manuscript) and S_n in terms of rho and omega
Delta2_simp:=(collect(algsubs(F_I=F_Ihat,algsubs(S_n=rho/omega,Delta2)),rho));
# Derivative of Delta2_simp wrt rho
dDelta2:=collect(simplify(diff(Delta2_simp,rho)),rho);
# only simplified numerator, I just want to see the quadratic
dDelta2_numer:=collect(numer(dDelta2)/gamma,rho);
# Note that I am assuming that the numerator of Delta' is in the form of (a*rho^2-b*rho+c)
a:=factor(coeff(dDelta2_numer,rho,2));
b:=-factor(coeff(dDelta2_numer,rho,1)); #Note the negative here!
# Maple doesn't accept c, so I put cc but in the manuscript it is shown by c
cc:=factor(coeff(dDelta2_numer,rho,0));
# Root finding of the quadratic: 
# Note that delta is always non-negative.
delta:=factor(simplify(b^2-4*a*cc));
root1:=simplify(b-sqrt(delta))/(2*a);
root2 :=simplify(b+sqrt(delta))/(2*a);
factor(numer(root1));
factor(numer(root2));
collect(factor(simplify(eval(dDelta2_numer,rho=omega))),[theta[w],theta[c]]);
# at rho=1 which is max the diff Delta wrt rho is
collect(simplify(a-b+cc),omega);
simplify(a-b+cc);
simplify(taylor(a-b+cc,omega=1,1));







# replacing rho with phi
F_Shat:=rho*omega/(omega-rho); # this is phi

Delta2_simp;
Delta3:=simplify(subs(rho=omega*phi/(omega+phi),Delta2_simp));
collect(numer(Delta3),phi);
factor(coeff(collect(numer(Delta3),phi),phi,2));
factor(coeff(collect(numer(Delta3),phi),phi,1));


# Derivative:
dDelta3:=collect(simplify(diff(Delta3,phi)),phi);



collect(simplify(diff(Delta3,phi)),phi);
dDelta3_numer:=collect(numer(collect(simplify(diff(Delta3,phi)),phi))/(gamma*omega),phi);
a3:=factor(coeff(dDelta3_numer,phi,2));
collect(a3,[theta[w],theta[c]]);
b3:=factor(coeff(dDelta3_numer,phi,1));
c3:=factor(coeff(dDelta3_numer,phi,0));
delta3:=factor(simplify(b3^2-4*a3*c3));

















