## create unified data frame
library(ggplot2)
library(dplyr)
library(tidyr)
library(directlabels)

source("SIRfunctions.R")
source("params.R")

#choose weight W_S 0.3 for TTI testing plan and 1 for the random testing.
W_S_random <- 1
W_S_targeted <- 0.3 

params_dat<-function(params){
  unpack(as.list(params))
  n_out <- 5 ## facets (rows/cols)
  n_in <- 41 ## grid N within facets
  tol <- 1e-10 ## to resolve the issue of very small numbers 
  # specify the ranges, FIXME, not to be hard coded!
  # eta_w <- seq(0,1, length.out=n_out)
  # eta_c <- seq(0,1, length.out=n_out)  
  rho <- seq(0,0.5, length.out=n_in) ## 0,0.01
  omega <- 0.5+tol #seq(0.5+tol,2 , length.out=n_in) #note omega must be > rho
  
  df1 <- expand.grid(N0=N0,beta=beta,gamma=gamma,omega=omega,rho=rho,
                     W_S=W_S,W_I=W_I,W_R=W_R,
                     p_S=p_S,p_I=p_I,p_R=p_R,
                     eta_w=eta_w,eta_c=eta_c)
  ## principle eigenvector
  eigvec <- t(apply(df1,1,function(params_in)eigvec_max(params=params_in)))  
dfout <-(df1 %>% 
            dplyr::mutate(
            R0=apply(df1,1,function(params_in)R0(params=params_in)),  
            R0_sub=ifelse(eta_w<eta_c, NA, R0), 
            theta_w=1-eta_w,
            theta_c=1-eta_c,
            Delta=ifelse(eta_w<eta_c, NA, 1-(R0*gamma/beta) ),
            Delta=ifelse(abs(Delta)<tol,0,Delta),
            I_u=eigvec[,1],
            I_n=eigvec[,2],
            I_p=eigvec[,3],
            I_c=eigvec[,4]
            ))
  return(dfout)
  }

df_random <- params_dat(params = update(params,W_S=W_S_random,eta_w=1,eta_c=0))
df_targeted <- params_dat(params = update(params,W_S=W_S_targeted,eta_w=1,eta_c=0))

out2 <- df_targeted %>%
  pivot_longer(c(I_u,I_n,I_p,I_c), names_to = "compartment", values_to = "value")
gg1 <- (ggplot(data=out2, aes(x=rho,y=value,col=compartment))
        + theme_bw()
        + geom_point()
        + ylab(label="Population fraction at DFE")
        + xlab(label="rho")
        )
direct.label(gg1,"last.bumpup")


##plot of R0 wrt rho, notice the increase
ggplot(data=out2, aes(x=rho,y=R0_sub))+geom_point()

