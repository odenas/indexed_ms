library(tidyverse)
library(knitr)

read_ds <- function(fname, set_input_type = FALSE){
  ds <- read_csv(fname)
  if(set_input_type){
    stopifnot("b_path" %in% colnames(ds))
    ds <- (ds %>% separate(b_path, into=c("inp_type", "sl", "tl", "alp", "mut_rate")))
  }
  ds
}


sorted_bpaths <- function(b_path){
  b_path <- as.character(levels(factor(b_path)))
  data.frame(b_path = b_path,
             mu_nr = as.numeric(simplify2array(lapply(strsplit(b_path, "_", fixed=TRUE), function(l) ifelse(length(l) == 4, l[[4]], NA)))),
             inp_type = sapply(b_path, function(s) ifelse(substr(s, 5, 6) == "1G", "G", "M")))
}


performance_boxplot <- function(p, tit, ylab_){
  p + geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha=0.1, width=0.1, size = 2) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
          legend.position="bottom") +
    labs(title = tit) + ylab(ylab_)
}


# input_description <- function(slen, tlen, itype, mperiod){
#   frm_len <- function(n){
#     if (n < 1000) sprintf("%d", n)
#     else if (n < 100000) sprintf("%dKB", round(n / 1000))
#     else if (n < 100000000) sprintf("%dMB", round(n / 1000000))
#     else sprintf("%dGB", round(n / 1000000000))
#   }
#
#   if(substr(itype, 1, 3) == "rnd")
#     sprintf("random with |s| = %s and |t| = %s", frm_len(slen), frm_len(tlen))
#   else
#     sprintf("mutated (every %d chars) with |s| = %s and |t| = %s", mperiod, frm_len(slen), frm_len(tlen))
# }
#



