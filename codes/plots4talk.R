# Making plots for presentations/slides etc in the "path"
path <- "~/projects/aliworkstation/references/presentations/epimodeling/pix/"
# R code is at: ~/projects/SIRmodel/codes/
# setwd("~/projects/SIRmodel/codes/")

source("SIRfunctions.R")
source("params.R")

# #################################
# Load the data
# #################################
df_random <- read.csv(file = 'random_test_df.csv') #random testing scenario
df_targeted <- read.csv(file = 'targeted_test_df.csv') #targeted testing scenario
#choose weight W_S 0.3 for TTI testing plan and 1 for the random testing.
W_S_random <- unique(df_random$W_S)
W_S_targeted <- unique(df_targeted$W_S)

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
mn <- min(df_random$R0_sub,df_targeted$R0_sub, na.rm = T)
mx <- max(df_random$R0_sub,df_targeted$R0_sub, na.rm = T)

# show_contours_1, width=24,height=24,message=FALSE,warning=FALSE}
brks_vec <- seq(0.8,1.05,by=0.05) # Break vector for unifying the legends in Random and TTI testing cases. 

# #################################
#  SubPlots the Random Testing Scenario:
# #################################
# Random Testing, Perfect isolation

# FIXME: df <- paste("df",scenario,sep="_")
# scenario <-"random"
# df <- df_random
scenario <-"targeted"
df <- df_targeted

# isolation for awaiting people
eta_w_val <-0
# eta_w_val <- unique(df_random$eta_w)

# isolation for confirmed people
eta_c_val <-0

df_temp <- df %>% dplyr::filter(eta_w %in% eta_w_val & eta_c %in% eta_c_val)
# no legend for Random plots
leg_pos <- "none" #ifelse(scenario=="targeted","right","none")
# Different title according to scenario
title_set <- ifelse(scenario=="targeted",
                    TeX(sprintf("Targeted Testing; W_S=%.1f, W_I=W_R=1",W_S_targeted )),
                    TeX(sprintf("Random Testing; W_S=W_I=W_R=1")))
  
p1 <- ggplot(df_temp,aes(x=omega,y=rho,z=R0_sub))+ theme_bw(base_size = 12) +
  xlab(TeX('$\\omega$, rate of test return  (1/day)')) +
  ylab(TeX('$\\rho$, testing intensity (1/day per capita)')) +
  theme(panel.spacing=grid::unit(0,"lines"),legend.position = leg_pos)
         

p1_temp <- (p1
            + geom_contour_filled(breaks=brks_vec)
            + geom_contour(breaks=1,alpha=0.5,colour="black")
            + facet_grid(eta_w~eta_c, labeller=label_special)
            + scale_x_continuous(expand=expansion(c(0,0)), n.breaks=3)
            + scale_y_continuous(expand=expansion(c(0,0)), n.breaks=3)
            + scale_fill_viridis_d(name=parse(text="R[0]"),drop=FALSE)
            + ggtitle(title_set)
            # + ggtitle(TeX(sprintf("Random Testing; W_S=W_I=W_R=1")))
            # + ggtitle(TeX(sprintf("Targetted Testing; W_S=%.1f, W_I=W_R=1",W_S_targeted )))
            + geom_rect(data=df_temp, fill=ifelse(df_temp$eta_c > df_temp$eta_w,"grey90","NA" ),
                        color= NA,
                        ymin=-1,
                        ymax=10,
                        xmin=-1,
                        xmax=2)
)
p1_temp

ggsave(p1_temp,
       filename = paste("R0_",scenario,"_","w",eta_w_val,"c",eta_c_val,".pdf",sep = "") ,
       width = 10, height = 10, units = "cm",
       path=path)

# Saving the column panels
ggsave(p1_temp,
       filename = paste("R0_",scenario,"_","w","all","c",eta_c_val,".pdf",sep = "") ,
       width = 12, height = 12, units = "cm",
       path=path)




