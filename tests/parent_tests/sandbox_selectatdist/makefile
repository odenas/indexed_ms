
EXEC =../../../fast_ms/bin/sandbox_selectatdist.x
DDIR = ../../datasets/big_paper2
R    = R CMD BATCH -q --vanilla --no-timing
#INPUT = $(DDIR)/rep_1000000s_dis_200000t_abcd_sim1000.s
#                rep_100000000s_dis_500000t_abcdefghijklmnopqrst_sim1000

INPUT  = $(notdir $(subst .s,,$(wildcard $(DDIR)/rep_100000000s_d*_sim1000.s) $(wildcard $(DDIR)/rnd_100000000s_d*_sim1000.s)))
TARGETS = $(addsuffix .csv,$(INPUT))
.PHONY : all rds

all : $(TARGETS)

rds : $(addsuffix _data.rds,all group)

group_data.rds : all_data.rds
	echo "library(tidyverse); saveRDS(readRDS('$<') %>% group_by(dist, method, inp_type, alp) %>% summarise(time_avg = mean(time_micro), time_sd = sd(time_micro)), '$@')" | $(R) /dev/stdin /dev/stdout

all_data.rds : $(TARGETS) save.R
	echo "OUT_PATH<-'$@';" | cat /dev/stdin save.R | $(R) /dev/stdin /dev/stdout


%.csv : $(DDIR)/%.s $(EXEC)
	$(EXEC) -load_maxrep 1 -load_cst 1  -s_path $< >$@

clean : 
	rm -f $(TARGETS)
