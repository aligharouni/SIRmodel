library(ggplot2)
library(dplyr)
library(tidyr)
library(directlabels)
library(latex2exp)
library(rootSolve)
library(McMasterPandemic)
# library(matlib) # used for matrix inverse function

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

# make dataframe 
df_random <- make_params_dat(params = params,
                eta_ws=0,eta_we=1, ## so theta_w
                eta_cs=0,eta_ce=1, ## so theta_c
                omega_s=0.1,omega_e=2,
                rho_s=0,rho_e=0.01)

df_targeted <- make_params_dat(params = update(params,W_S=W_S_targeted),
                eta_ws=0,eta_we=1, ## so theta_w
                eta_cs=0,eta_ce=1, ## so theta_c
                omega_s=0.1,omega_e=2,
                rho_s=0,rho_e=0.01)
## high rho and omega
df_random_h <- make_params_dat(params = params,
                             eta_ws=0,eta_we=1, ## so theta_w
                             eta_cs=0,eta_ce=1, ## so theta_c
                             omega_s=0.25,omega_e=2,
                             rho_s=0,rho_e=0.25-0.000001)

df_targeted_h <- make_params_dat(params = update(params,W_S=W_S_targeted),
                               eta_ws=0,eta_we=1, ## so theta_w
                               eta_cs=0,eta_ce=1, ## so theta_c
                               omega_s=0.25,omega_e=2,
                               rho_s=0,rho_e=0.25-0.000001)

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

n_contour <- 8 ## number of desired contours or bins
bins <- cut(sort(unique(c(df_random$Delta,df_targeted$Delta))), 
            n_contour)
brks <- levels(bins)
brks_vec <- unique(as.numeric(unlist(lapply(strsplit(brks, ","), function(x) gsub("\\(|]", "", x)))))
brks_vec <- sort(ifelse(brks_vec[]< 0,0,brks_vec[]), decreasing = TRUE) ## replace neg with 0 and reorder

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
            # + geom_contour_filled(breaks=brks)
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
# 3. Plots for high rho Scenario:
# #################################
bins <- cut(sort(unique(c(df_random_h$Delta,df_targeted_h$Delta))), 
            n_contour)
brks <- levels(bins)
brks_vec <- unique(as.numeric(unlist(lapply(strsplit(brks, ","), function(x) gsub("\\(|]", "", x)))))
brks_vec <- sort(ifelse(brks_vec[]< 0,0,brks_vec[]), decreasing = TRUE) ## replace neg with 0 and reorder

df_temp2 <- df_random_h

p12 <- (ggplot(df_temp2,aes(x=1/omega,y=rho,z=Delta))
       + theme_bw()
       + xlab(TeX('$\\1/omega$, mean test return time (day)'))
       + ylab(TeX('$\\rho$, testing intensity (1/day per capita)'))
)
## There might be a solution for the tick label collision with expand_limits(), but that's hard to do across facets
## see https://stackoverflow.com/questions/41575045/avoiding-axis-tick-label-collision-in-faceted-ggplots
p1_temp2 <- (p12
            # + geom_contour_filled(breaks=brks)
            + geom_contour_filled(breaks=brks_vec)
            + geom_contour(breaks=threshold,alpha=0.5,colour="black")
            + facet_grid(theta_w~theta_c, labeller=label_special) 
            + scale_x_continuous(expand=expansion(c(0,0)), n.breaks=3)
            + scale_y_continuous(expand=expansion(c(0,0)), n.breaks=3)
            + scale_fill_viridis_d(name=parse(text="Delta"),drop=FALSE)
            + geom_rect(data=df_temp2, fill=ifelse(df_temp2$theta_c < df_temp2$theta_w,"grey90","NA" ),
                        color= NA,
                        ymin=-1,
                        ymax=10,
                        xmin=-1,
                        xmax=10)
            + theme(panel.spacing=grid::unit(0,"lines"))
)
p1_temp2
# #################################
# 4. Plot the Random Testing Scenario:
# #################################

ggsave(p1_temp2 + ggtitle(TeX("w_S=w_I=w_R=1")) +
         theme(legend.position = "none"),
       filename = "R0contour_random2.pdf" ,
       width = 12, height = 12, units = "cm")

# #################################
# 5. Plot the Targeted Testing Scenario:
# #################################

ggsave((p1_temp2 %+% df_targeted_h) +
         ggtitle(TeX(sprintf("w_S=%.1f, w_I=w_R=1",W_S_targeted))) +
         theme(legend.position = c(0.2, 0.3),legend.text = element_text(size = 8)),
       filename = "R0contour_TTI2.pdf" ,
       width = 12, height = 12, units = "cm")

