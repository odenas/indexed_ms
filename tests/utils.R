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

PEAK_SPACE_LAB <- c("alg", "bwt_bwt", "bwt_alp", "bwt_wtree",
                    "stree_bpsupp", "stree_bpsupp_select", "stree_bpsupp_rank", "ms", "runs")
SYS_SPACE_LAB <- c("resident_mem", "virtual_mem")


size_s_f <- function(x)
  as.numeric(simplify2array(strsplit(x, "_"))[2,])

size_t_f <- function(x)
  as.numeric(sapply(simplify2array(strsplit(x, "_"))[1,], function(ss) substring(ss, 14)))

peak_space <- function(ds_, label){
  item_set <- c("alg", "bwt_bwt", "bwt_alp", "bwt_wtree",
                "stree_bpsupp", "stree_bpsupp_select", "stree_bpsupp_rank",
                "ms", "runs")
  (filter(ds_, item %in% item_set, measuring=="space") %>%
      mutate(label = label, size_s = size_s_f(prefix)))
}

system_space <- function(ds_, label)
  (filter(ds_, item %in% c("resident_mem", "virtual_mem"), measuring=="space") %>%
     mutate(label = label, size_s = size_s_f(prefix)))

peak_time <- function(ds_, label)
  (filter(ds_, measuring=="time", item=="total_time") %>%
     mutate(label = label, size_s = size_s_f(prefix), size_t = size_t_f(prefix)))



linear_fit_plot <- function(ds, fit, xlab, ylab){
  ggplot(ds, aes(size_s, value)) + geom_point() + geom_line() +
    labs(title=sprintf("%.2f + %.4f * |s|",
                       coef(fit)[1], coef(fit)[2]), x=xlab, y=ylab) +
    geom_abline(intercept=coef(fit)[1], slope=coef(fit)[2], color='blue', alpha=0.5) +
    geom_abline(intercept=0, slope=1, color='green', alpha=0.5)
}

bwt_usage <- function(ds_){
  filter(ds_, item %in% c("bwt_wtree", "bwt_bwt", "bwt_alp")) %>%
    group_by(prefix, label) %>% summarize(size_s=max(size_s), value=sum(value))
}

stree_usage <- function(ds_){
  st_size_items <- c("stree_bp", "stree_bpsupp",
                     "stree_bpsupp_select", "stree_bpsupp_rank")
  filter(ds_, item %in% st_size_items) %>%
    group_by(prefix, label) %>%
    summarize(size_s=max(size_s), value=sum(value))
}

vec_usage <- function(ds_){
  filter(ds_, item %in% c("runs", "ms")) %>%
    group_by(prefix, label) %>% summarize(size_s=max(size_s), value=sum(value))
}

