library(tidyverse)
load_dsets <- function(dset_dir, cache, refresh=TRUE){
  if((!file.exists(cache)) || refresh){
    all_dt <- do.call(rbind, lapply(list.files(dset_dir, pattern = "*.csv"), function(fname){
      parts <- strsplit(fname, "_", fixed=TRUE)[[1]]
      
      read_csv(sprintf("%s/%s", dset_dir, fname)) %>% mutate(inp_type = sprintf("%s_%s", parts[1], parts[3]), alp=parts[5])
    }))
    message(sprintf("* dumping to %s", cache))
    saveRDS(all_dt, cache)
  }
  message(sprintf("* loading from %s", cache))
  readRDS(cache)
}
load_dsets(".", OUT_PATH, TRUE)
