#!/usr/bin/R Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least 1 arguments must be supplied (input file)", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = 0.01
  args[3] = "corr_plot"
}

if (!'data.table' %in% installed.packages()) install.packages("data.table") 
  library(data.table)
if (!'ggplot2' %in% installed.packages()) install.packages("ggplot2") 
  library(ggplot2)
if (!'ggpmisc' %in% installed.packages()) install.packages("ggpmisc") 
  library(ggpmisc)

# setwd('/home/rskick/Projects/y/')

plotting <- function(name, df, a){
  ggplot(df, aes(V1, V2))+
    geom_point(alpha=a)+
    stat_smooth(method = 'lm', formula = y ~ x, col='red', )+
    stat_poly_eq(formula = y ~ x, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~~~~")), 
                 parse = TRUE) +
    scale_x_continuous(name = 'AF')+
    scale_y_continuous(name = 'GNOMAD AF')+
    theme_bw()
  ggsave(sprintf("%s.pdf", name))
  }

df <- fread(args[1])
df0_01 <- df[V2 > 0.01]

plotting(sprintf("%s_all", args[3]), df, args[2])
plotting(sprintf("%s_freq", args[3]), df0_01, args[2])

  