## create unified data frame
library(ggplot2)
library(dplyr)
library(tidyr)
library(directlabels)

source("SIRfunctions.R")
source("params.R")

sup_df_random <- make_params_dat(params = params,
                             eta_ws=1,eta_we=1, ## so theta_w=0
                             eta_cs=0,eta_ce=0, ## so theta_c=1
                             omega_s=0.5001,omega_e=0.5001,
                             rho_s=0,rho_e=0.5)

sup_df_targeted <- make_params_dat(params = update(params,W_S=0.3),
                             eta_ws=1,eta_we=1, ## so theta_w=0
                             eta_cs=0,eta_ce=0, ## so theta_c=1
                             omega_s=0.5001,omega_e=0.5001,
                             rho_s=0,rho_e=0.5)
out2 <- sup_df_targeted %>%
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

