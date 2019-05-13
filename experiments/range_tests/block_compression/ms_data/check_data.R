library(tidyverse)
library(tidylog)
rm(list=ls())

L <- 6264717366
dd <- (
  read_csv("all.csv")
  %>% separate(compr, c("compr1", "compr"))
  %>% mutate(compr=if_else(is.na(compr), "none", compr),
             block_size=block_end - block_start,
             full = if_else(block_size == L, "full", "block"))
  %>% select(-compr1)
)
dd
full_d <- (dd
           %>% filter(block_size == 6264717366)
           %>% select(-block_start, -block_end)
           %>% mutate(compr_rate = compr_size / block_size))
full_d

# how does compression rate change along the ms vector?
ggplot(dd, aes(block_start, compr_size/block_size, color=compr, shape=full)) +
  geom_point() + geom_line()

dd <- dd %>% inner_join(full_d %>% rename(tot_byte_size=compr_size, tot_len=block_size), by=c("compr"))
dd

(dd %>% group_by(compr, full)
  %>% summarise(compr_size = sum(compr_size))
  %>% spread(full, compr_size) %>% mutate(compr_r = block / full)
)
