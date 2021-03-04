library(shellpipes)
rpcall("params.Rout params.R SIRfunctions.rda")

loadEnvironments()

# Parameters used in SIR model
# Path: ~/projects/SIRmodel/codes

params <- c(
N0=1000000, #total population size
beta=.8, #transmission rate
gamma=1/3, #recovery rate
omega=0.25, #test returning rate 
rho=0.01, #testing rate percapita (also used 0.8, 1/3)
W_S=1, W_I=1, W_R=1, #Testing relative weights
p_S=0, p_I=1, p_R=0.5, #test specificity, i.e., prob of waiting for being tested positive
eta_w=0.02, eta_c=0.01, #isolation parameter (eta=0 is the perfect isolation) 
s=2 # isolation ratio parameter, i.e., s=eta_w/eta_c
            )
            
class(params) <- "params_pansim"

state_init <- c(S_u=params[["N0"]], S_n=0,
                I_u=1,I_n=0,I_p=0,I_t=0,
                R_u=0,R_n=0,R_p=0,R_t=0,
                N=0,P=0)
                


state_dfe <- c(S_u=Su_dfe(params), S_n=Sn_dfe(params),
                I_u=0,I_n=0,I_p=0,I_t=0,
                R_u=0,R_n=0,R_p=0,R_t=0,
                N=0,P=0)

saveEnvironment()
