library(ggplot2)
library(dplyr)
library(tidyr)
library(directlabels)
library(latex2exp)
library(rootSolve)
library(McMasterPandemic)

library(shellpipes)

rpcall("sir_plot.Rout sir_plot.R params.rda")

loadEnvironments()

## setwd("/home/ag/projects/SIRmodel/codes/")

# #################################
# Making Dataframe part
# #################################
# R0 contour plots
n_out <- 5 ## facets (rows/cols)
n_in <- 41 ## grid N within facets

#choose weight W_S 0.3 for TTI testing plan and 1 for the random testing.
W_S_random <- 1
W_S_targeted <- 0.3 

# make dataframe and save it, if it doesn't exists already.
df_random <- eval_R0(params = update(params,W_S=W_S_random))
df_targeted <- eval_R0(params = update(params,W_S=W_S_targeted))

##important contour, ie R0=1 thus threshold=1 when plotting R0 contours, or Delta(R0=1)
threshold <- 1-(params[["gamma"]]/params[["beta"]]) ## Or 0?  corresponding to R0=1 
# #################################
# Plotting Part
# #################################
# Setting the plots
## hacked (not very robustly) to chain label_both() and label_parsed() ...
label_special <- function (labels, multi_line=FALSE, sep = "== ") {
  value <- ggplot2:::label_value(labels, multi_line = multi_line)
  variable <- ggplot2:::label_variable(labels, multi_line = multi_line)
  variable[[1]] <- gsub("_(.)$","[\\1]",variable[[1]])
  ## not using multiple faceting variables on each margin, so forget paste
  out <- Map(paste, variable, value, sep = sep)
  out[[1]] <- lapply(out[[1]], function(expr) c(parse(text=expr)))
  out
}

# verify the break range in brks_vec (FIXME)
mn <- min(df_random$Delta,df_targeted$Delta, na.rm = T)
mx <- max(df_random$Delta,df_targeted$Delta, na.rm = T) ## take mn and mx and set the brks_vec, smarter way?

# show_contours_1, width=24,height=24,message=FALSE,warning=FALSE}
brks_vec <- seq(.165,0,by=-0.03) # Break vector for unifying the legends in Random and TTI testing cases. 

df_temp <- df_random
# df_temp <- df_targeted
p1 <- (ggplot(df_temp,aes(x=1/omega,y=rho,z=Delta))
    + theme_bw()
    + xlab(TeX('$\\1/omega$, mean test return time (day)'))
    + ylab(TeX('$\\rho$, testing intensity (1/day per capita)'))
)
## There might be a solution for the tick label collision with expand_limits(), but that's hard to do across facets
## see https://stackoverflow.com/questions/41575045/avoiding-axis-tick-label-collision-in-faceted-ggplots
p1_temp <- (p1
            + geom_contour_filled(breaks=brks_vec)
            + geom_contour(breaks=threshold,alpha=0.5,colour="black")
            + facet_grid(theta_w~theta_c, labeller=label_special) 
            + scale_x_continuous(expand=expansion(c(0,0)), n.breaks=3)
            + scale_y_continuous(expand=expansion(c(0,0)), n.breaks=3)
            + scale_fill_viridis_d(name=parse(text="Delta"),drop=FALSE)
            + geom_rect(data=df_temp, fill=ifelse(df_temp$theta_c < df_temp$theta_w,"grey90","NA" ),
                        color= NA,
                        ymin=-1,
                        ymax=10,
                        xmin=-1,
                        xmax=10)
    + theme(panel.spacing=grid::unit(0,"lines"))
)
p1_temp
# #################################
# 1. Plot the Random Testing Scenario:
# #################################

ggsave(p1_temp + ggtitle(TeX("w_S=w_I=w_R=1")) +
       theme(legend.position = "none"),
       filename = "R0contour_random.pdf" ,
       width = 12, height = 12, units = "cm")

# #################################
# 2. Plot the Targeted Testing Scenario:
# #################################

ggsave((p1_temp %+% df_targeted) +
       ggtitle(TeX(sprintf("w_S=%.1f, w_I=w_R=1",W_S_targeted))) +
       theme(legend.position = c(0.2, 0.3),legend.text = element_text(size = 8)),
       filename = "R0contour_TTI.pdf" ,
       width = 12, height = 12, units = "cm")

# #################################
# 3. Plot R0 as a function of rho (testing intensity):
# #################################
eta_w <- 0.5
eta_c <- eta_w
gamma <- 1/4
W_S <- 1
omega <- 10^6
rho <- seq(0.5,1, length.out=100)
f <- function(params){
  unpack(as.list(params))
  return(gamma*(W_S/W_I)*1/(sqrt(eta_w/gamma)-1))
}
f(params = update(params,c(eta_w=eta_w,eta_c=eta_c,
                           gamma=gamma)))

df1 <- expand.grid(N0=params[["N0"]],beta=params[["beta"]],gamma=gamma,
                   omega=omega,rho=rho,
                   W_S=W_S,W_I=params[["W_I"]],W_R=params[["W_R"]],
                   p_S=params[["p_S"]],p_I=params[["p_I"]],p_R=params[["p_R"]],
                   eta_w=eta_w,eta_c=eta_c)
df2<- data.frame(df1, R0=apply(df1,1,function(params_in)R0(params=params_in)),
                      f=apply(df1,1,function(params_in)f(params = params_in)),
                      S_u=apply(df1,1,function(params_in)Su_dfe(params = params_in))
                )
ggplot(df2,aes(x=rho,y=R0))+geom_point()
