library(shellpipes)

# Functions for the SIR model
# By: A Gharouni

# install remotes package if necessary:
while (!require(remotes)) {
  install.packages("remotes")
}
## install development version of bbmle:
if (!require("bbmle") || packageVersion("bbmle") < "1.0.23.5") {
  remotes::install_github("bbolker/bbmle")
}
## install the target package and all its dependencies:
while (!require(McMasterPandemic)) {
  remotes::install_github("bbolker/McMasterPandemic",
                          dependencies = TRUE,
                          build_vignettes = TRUE
  )
}

unpack <- McMasterPandemic::unpack

# Testing functions
sigma <- function(state,params){
  unpack(as.list(c(state,params)))
  W <- W_S*S_u+W_I*I_u+W_R*R_u
  testing_rate <- rho*N0/W
  return(testing_rate)
}

# SIR model
sir.model <- function(time,state,params){
  unpack(as.list(c(state,params)))
  ## Force of Infection  
  Lambda <- beta * (I_u + eta_w*(I_n+I_p) + eta_c*I_t)/N0
  ## scaling the weights
  # W <- W_S*S_u+W_I*I_u+W_R*R_u 
  # sigma <- rho*N0/W
  sigma <- sigma(state,params)
  
  # testing intensity
  F_S <- sigma*W_S
  F_I <- sigma*W_I
  F_R <- sigma*W_R
  # Equations
  dS_u.dt <- -Lambda*S_u - (1-p_S) * F_S * S_u + omega * S_n
  dS_n.dt <- -Lambda*S_n + (1-p_S) * F_S * S_u - omega * S_n 
  dI_u.dt <-  Lambda*S_u + omega * I_n - F_I * I_u - gamma * I_u
  dI_n.dt <-  Lambda*S_n + (1-p_I) * F_I * I_u  - omega * I_n - gamma * I_n
  dI_p.dt <-  p_I * F_I * I_u - omega * I_p - gamma * I_p 
  dI_t.dt <-  omega * I_p - gamma * I_t 
  dR_u.dt <- gamma * I_u + omega * R_n - F_R * R_u  
  dR_n.dt <- gamma * I_n + (1-p_R) * F_R * R_u - omega * R_n
  dR_p.dt <- gamma * I_p + p_R * F_R * R_u - omega * R_p
  dR_t.dt <- gamma * I_t + omega * R_p
  dN.dt <- omega * (S_n + I_n + R_n)
  dP.dt <- omega *(I_p + R_p)
  
  # return the rate of change
  dxdt <- c(dS_u.dt,dS_n.dt,dI_u.dt,dI_n.dt,dI_p.dt,dI_t.dt,dR_u.dt,dR_n.dt,dR_p.dt,dR_t.dt,dN.dt,dP.dt)
  ## }
  return(list(dxdt))
}

# TODO: fix the sigma the issue is when DFE(state,params), the multiroot() gives error.
DFE <- function(S,params){
  unpack(as.list(params))
  # S: Susceptibles
  S_u <- S[1]
  S_n <- S[2]
  I_u <-0
  R_u <-0
    #Weighted untested people
  W <- W_S*S_u+W_I*I_u+W_R*R_u
  sigma <- rho*N0/W
  F_S <- sigma*W_S
  
  return(c(eq1=S_u+S_n-N0,
           eq2=-F_S*S_u+omega*S_n))
}

# awaiting for negative tested Susceptible at DFE
Sn_dfe <-function(params){
  unpack(as.list(params))
  return((rho*(1-p_S)*N0)/omega)
}

# untested Susceptible at DFE
Su_dfe <-function(params){
  unpack(as.list(params))
  return(N0-Sn_dfe(params))
}

# Model simulation function
run.sir <- function(model, params,state,sim_time){
  # use update(params,beta=2)
  library(deSolve)
  unpack(as.list(c(state,params)))
  out <- as.data.frame(
    ode(
    func=model,
    y=state,
    times= sim_time, 
    parms=params
  ))
  return(out)
}


F_I <- function(state,params){
  # Returns F_I at the specified state
  unpack(as.list(c(state,params)))
  
  # sigma <- rho*N0/(W_S*Su) #at DFE 
  sig <- sigma(state,params)
  return(sig*W_I)
}

# Basic Reproduction Number
R0<-function(state=state_dfe,params){
  unpack(as.list(c(state,params)))
  Sn <- Sn_dfe(params)
  Su <- Su_dfe(params)
  Fi <- F_I(state,params) #at DFE
  
  A <- gamma*(omega+gamma+eta_w*Fi) + omega*eta_c*p_I*Fi
  B <- (omega+eta_w*(Fi+gamma))*gamma+(eta_w*gamma+eta_c*omega)*(omega*p_I*Fi)/(omega+gamma)
  C <- beta/(N0*gamma*(gamma*(omega+gamma)+Fi*(gamma+omega*p_I)))
  return((A*Su+B*Sn)*C)
}

# R0 expression with eta_w/eta_c factored
R02 <- function(state=state_dfe,params){
  unpack(as.list(c(state,params)))
  Sn <- Sn_dfe(params)
  Su <- Su_dfe(params)
  Fi <- F_I(state,params) #at DFE
  
  A2 <- gamma*(p_I*Fi*omega*Sn+(omega+gamma)*(gamma*Sn+Fi*N0))
  B2 <- p_I*Fi*omega*(gamma*Su+omega*N0)
  D <- gamma*(omega+gamma)*(gamma*Su+omega*N0)
  C2 <- beta/((N0*gamma*(gamma*(omega+gamma)+Fi*(gamma+omega*p_I)))*(omega+gamma))
  # s <- eta_w/eta_c
  return((eta_c*(s*A2+B2)+D)*C2)
}

# make grid and calc R0 for plotting
eval_R0 <- function(state=state_dfe,params){
  # input the params and their range, this function makes a grid dataframe, calls R0 and outputs the csv file
  unpack(as.list(c(state,params)))
  tol <- 1e-10 ## to resolve the issue of very small numbers 
  # specify the ranges, FIXME, not to be hard coded!
  beta  <- params[["beta"]]  ## set beta high enough to allow R0>1 in worst case
  eta_w <- seq(0,1, length.out=n_out)
  eta_c <- seq(0,1, length.out=n_out) #0.001 to 0.5 
  rho <- seq(0,0.01, length.out=n_in)
  omega <- seq(0.1,2 , length.out=n_in) #note omega must be > rho, was 0.011
  
  df1 <- expand.grid(N0=params[["N0"]],beta=beta,gamma=params[["gamma"]],omega=omega,rho=rho,
                     W_S=W_S,W_I=params[["W_I"]],W_R=params[["W_R"]],
                     p_S=params[["p_S"]],p_I=params[["p_I"]],p_R=params[["p_R"]],
                     eta_w=eta_w,eta_c=eta_c)
  df2<- data.frame(df1,
                   R0=apply(df1,1,function(params_in)R0(params=params_in)),
                   Fi= apply(df1,1,function(params_in)F_I(state=state, params=params_in)))
  # This is for plotting purposes:
  df2 <- (df2 %>% 
            dplyr::mutate(R0_sub=ifelse(eta_w<eta_c, NA, R0), 
                          Eta_w=1-eta_w,
                          Eta_c=1-eta_c,
                          Delta=ifelse(eta_w<eta_c, NA, 1-(R0*gamma/beta) ),
                          Delta=ifelse(abs(Delta)<tol,0,Delta)
                          ) 
          )
  return(df2)
}



# ###################
# Extra Stuff unused
# ###################

# Conditional stop for ODE desolver
# rootfun <- function (time, state,params) { 
#   unpack(as.list(c(state,params)))
#   return(I_u - 1) }

# ode(
#   func=sir.model,
#   y=state_init,
#   times= d, 
#   parms=params,
#   atol = 1e+1, rtol = 1e+1
#   rootfun = rootfun,
#   method="lsodar"
# )

saveEnvironment()

