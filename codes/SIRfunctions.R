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
remotes::install_github("bbolker/McMasterPandemic",
                        dependencies = TRUE,
                        build_vignettes = TRUE
)


# SIR model
sir.model <- function(time,state,params){
  unpack(as.list(c(state,params)))
  ## Force of Infection  
  Lambda <- beta * (I_u + eta_w*(I_n+I_p) + eta_c*I_t)/N0
  ## scaling the weights
  W <- W_S*S_u+W_I*I_u+W_R*R_u
  sigma <- rho*N0/W
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


DFE <- function(S,params){
  # S: Susceptibles
  S_u <- S[1]
  S_n <- S[2]
  I_u <-0
  R_u <-0
  unpack(as.list(params))
  W <- W_S*S_u+W_I*I_u+W_R*R_u #Weighted untested people
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




# ###################
# Extra Stuff unused
# ###################

# Conditional stop for desolver
# rootfun <- function (time, state,params) { 
#   unpack(as.list(c(state,params)))
#   return(I_u - 1) }

