rm(list=ls())

library(ggplot2)
library(dplyr)
library(knitr)
library(reshape2)

read_ds <- function(fname, label=NA){
  if(!is.na(label))
    read.csv(fname, header = TRUE, stringsAsFactors = FALSE) %>% mutate(label = label)
  else
    read.csv(fname, header = TRUE, stringsAsFactors = FALSE)
}

PEAK_SPACE_LAB <- c("alg",
                    "bwt_bwt", "bwt_alp", "bwt_wtree",
                    "stree_bpsupp", "stree_bpsupp_select", "stree_bpsupp_rank",
                    "stree_csa", "stree_bp", "stree_bpsupp",
                    "ms", "runs")
VEC_SPACE_LAB <- c("ms", "runs")
SYS_SPACE_LAB <- c("resident_mem", "virtual_mem")

linear_fit_plot <- function(ds, fit, xlab, ylab){
  ggplot(ds, aes(len_s, value)) + geom_point() + geom_line() +
    labs(title=sprintf("%.2f + %.4f * |s|",
                       coef(fit)[1], coef(fit)[2]), x=xlab, y=ylab) +
    geom_abline(intercept=coef(fit)[1], slope=coef(fit)[2], color='blue', alpha=0.5) +
    geom_abline(intercept=0, slope=1, color='green', alpha=0.5)
}

