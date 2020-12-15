
source("SIRfunctions.R")
source("params.R")

# #################################
# Making Dataframe part
# #################################
# R0 contour plots
n_out <- 5 ## facets (rows/cols)
n_in <- 41 ## grid N within facets

#choose weight W_S 0.3 for TTI testing plan and 1 for the random testing.
W_S_random <- 1
W_S_targeted <- 0.3 

# make dataframe and save it. 
# FIXME: I feel that this part needs to be separated from here?
# eval_R0(params = update(params,W_S=W_S_random),filename="random_test_df.csv")
# eval_R0(params = update(params,W_S=W_S_targeted),filename="targeted_test_df.csv")

# #################################
# Load the data
# #################################
df_random <- read.csv(file = 'random_test_df.csv')
df_targeted <- read.csv(file = 'targeted_test_df.csv')

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

df_temp <- df_random
p1 <- ggplot(df_temp,aes(x=omega,y=rho,z=R0_sub))+ theme_bw() +
  xlab(TeX('$\\omega$, Returning test result rate (1/day)')) +
  ylab(TeX('$\\rho$, testing intensity (1/day per capita)')) +
  # theme(panel.spacing=grid::unit(0,"lines"),legend.position = c(0.89, 0.73)) # initial theme, comment it for TTI testing plot
  theme(panel.spacing=grid::unit(0,"lines"),legend.position = "none") #AG, uncomment for Random Case
p1_temp <- (p1
            + geom_contour_filled(breaks=brks_vec)
            + geom_contour(breaks=1,alpha=0.5,colour="black")
            + facet_grid(eta_w~eta_c, labeller=label_special)
            + scale_x_continuous(expand=expansion(c(0,0)), n.breaks=3)
            + scale_y_continuous(expand=expansion(c(0,0)), n.breaks=3)
            + scale_fill_viridis_d(name=parse(text="R[0]"),drop=FALSE)
            + geom_rect(data=df_temp, fill=ifelse(df_temp$eta_c > df_temp$eta_w,"grey90","NA" ),
                        color= NA,
                        ymin=-1,
                        ymax=10,
                        xmin=-1,
                        xmax=2)
)

# #################################
# 1. Plot the Random Testing Scenario:
# #################################

ggsave(p1_temp + ggtitle(TeX("W_S=W_I=W_R=1")),
       filename = "R0contour_random.pdf" ,
       width = 14, height = 14, units = "cm",
       path="../pix/")

# #################################
# 2. Plot the Targeted Testing Scenario:
# #################################

ggsave((p1_temp %+% df_targeted) +
        ggtitle(TeX(sprintf("W_S=%.1f, W_I=W_R=1",W_S_targeted ))),
    filename = "R0contour_TTI.pdf" ,
    width = 14, height = 14, units = "cm",
    path="../pix/")





