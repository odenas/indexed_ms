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

sorted_bpaths <- function(b_path){
  b_path <- as.character(levels(factor(b_path)))
  data.frame(b_path = b_path,
             mu_nr = as.numeric(simplify2array(lapply(strsplit(b_path, "_", fixed=TRUE), function(l) ifelse(length(l) == 4, l[[4]], NA)))),
             inp_type = sapply(b_path, function(s) ifelse(substr(s, 5, 6) == "1G", "G", "M")))
}

read_time_ds <- function(path = 'lazy_vs_nonlazy_data/lazy_vs_nonlazy.csv',
                         items=c("ms_bvector", "runs_bvector")){
  (read_ds(path) %>%
     filter(measuring == "time", item %in% items) %>%
     select(len_s, len_t, value, item, label, b_path))
}

read_lazy_call_cnt <- function(path = 'lazy_vs_nonlazy_data/lazy.csv'){
  lazy_calls <- sapply(0:3000, function(i) sprintf("consecutive_lazy_wl_calls%d", i))
  ds <- (read_ds(path) %>%
           filter(item %in% lazy_calls, label=="lazy") %>%
           mutate(count = value) %>%
           group_by(b_path, item) %>%
           summarise(len_s = mean(len_s),
                     len_t = mean(len_t),
                     value = mean(value),
                     count = mean(count)))
  ds$nwlcalls = as.numeric(sapply(ds$item, function(s) substr(s, 26, 300)))
  ds
}

