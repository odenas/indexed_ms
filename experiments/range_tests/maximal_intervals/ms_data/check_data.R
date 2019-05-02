library(tidyverse)
library(tidylog)
rm(list=ls())

dd <- read_csv('runs-1000.csv')
dd
ggplot(dd, aes(start, cnt)) + geom_point(alpha=0.4)
ggplot(dd, aes(start, (cnt0 - cnt1) / cnt)) + geom_point()

ggplot(dd, aes(factor(start), idx)) + geom_boxplot()
(dd %>% group_by(start) %>% summarise(n=n()))
ggplot(dd %>% group_by(start) %>% summarise(len=mean(len), idx=mean(idx) ), aes(start, idx)) + geom_line()

table(
  dd
  %>% group_by(start)
  %>% mutate(a = c(0, value[1:n() - 1] - value[2:n()]))
  %>% mutate(a = (a!= 0))
  %>% select(start, a)
)

table(dd %>% select(-idx) %>% group_by(start) %>% ungroup())


ggplot(dd, aes(idx, value, color=factor(start))) + geom_line()
