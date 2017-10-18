
EXEC =../../../linux_build/sandbox_test_maxrep/main
DDIR = ../../datasets/small_paper3
R    = R CMD BATCH -q --vanilla --no-timing

INPUT  = $(notdir $(subst .s,,$(wildcard $(DDIR)/rep_1000000s_dis_*.s) $(wildcard $(DDIR)/rnd_1000000s_dis*.s)))
TARGETS = $(addsuffix .csv,$(INPUT))
.PHONY : all rds

all : $(TARGETS)

rds : $(addsuffix .rds,maximal0_close0 maximal0_close1 maximal1_close0 maximal1_close1)

maximal0_close0.rds : all_data.rds
	echo "library(tidyverse); saveRDS(readRDS('$<') %>% filter(is_maximal == 0, close == 0), '$@')" | $(R) /dev/stdin /dev/stdout

maximal0_close1.rds : all_data.rds
	echo "library(tidyverse); saveRDS(readRDS('$<') %>% filter(is_maximal == 0, close == 1), '$@')" | $(R) /dev/stdin /dev/stdout

maximal1_close0.rds : all_data.rds
	echo "library(tidyverse); saveRDS(readRDS('$<') %>% filter(is_maximal == 1, close == 0), '$@')" | $(R) /dev/stdin /dev/stdout

maximal1_close1.rds : all_data.rds
	echo "library(tidyverse); saveRDS(readRDS('$<') %>% filter(is_maximal == 1, close == 1), '$@')" | $(R) /dev/stdin /dev/stdout

all_data.rds : $(TARGETS) save.R
	echo "OUT_PATH<-'$@';" | cat /dev/stdin save.R | $(R) /dev/stdin /dev/stdout

%.csv : $(DDIR)/%.s $(EXEC)
	$(EXEC) -load_maxrep 1 -load_cst 1  -s_path $< -t_path $(subst .s,.s,$<) >$@

clean : 
	rm -f $(TARGETS)
