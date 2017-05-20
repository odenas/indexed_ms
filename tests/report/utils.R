rm(list=ls())

library(tidyverse)
library(knitr)

read_ds <- function(fname, set_input_type = FALSE){
  parse_b_path <- function(x){
    if(length(x) == 4) sprintf("mut_%s", x[[4]]) else "rnd"
  }

  ds <- read_csv(fname)
  if(set_input_type){
    ds$inp_type <- simplify2array(lapply(strsplit(ds$b_path, "_", fixed = TRUE), parse_b_path))
  }
  ds
}


sorted_bpaths <- function(b_path){
  b_path <- as.character(levels(factor(b_path)))
  data.frame(b_path = b_path,
             mu_nr = as.numeric(simplify2array(lapply(strsplit(b_path, "_", fixed=TRUE), function(l) ifelse(length(l) == 4, l[[4]], NA)))),
             inp_type = sapply(b_path, function(s) ifelse(substr(s, 5, 6) == "1G", "G", "M")))
}

read_lazy_call_cnt <- function(path = '../lazy_vs_nonlazy_data/lazy.csv'){
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




