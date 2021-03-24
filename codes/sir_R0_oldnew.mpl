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
C_:= 1/((omega+gamma)*(gamma*omega+gamma^2+F_I*gamma+omega*p_I*F_I));
A_ * theta[w] + B_ * theta[c];
Delta:=%*C_;
R0_:= simplify( algsubs(S_u+S_n=1,beta/gamma * (1-Delta)) );
simplify(algsubs(theta[c]=1-eta[c],algsubs(theta[w]=1-eta[w],R0_)));
# Check between old and new notation; checked!
simplify(%-subs(N_0=1,algsubs(S_u+S_n=N_0,R0)));





