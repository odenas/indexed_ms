
EXEC =../../../linux_build/sandbox_test_maxrep/main
DDIR = ../../datasets/big_paper2
#INPUT = $(DDIR)/rep_1000000s_dis_200000t_abcd_sim1000.s
#                rep_100000000s_dis_500000t_abcd_sim1000

INPUT  = $(notdir $(subst .s,,$(wildcard $(DDIR)/rep_100000000s_dis*_sim1000.s) $(wildcard $(DDIR)/rnd_100000000s_dis*_sim1000.s)))
TARGETS = $(addsuffix .csv,$(INPUT))

all : all_data.rds

all_data.rds : $(TARGETS) save.R
	echo "OUT_PATH<-'$@';" | cat /dev/stdin save.R | R CMD BATCH -q --vanilla --no-timing /dev/stdin /dev/stdout

%.csv : $(DDIR)/%.s $(EXEC)
	$(EXEC) -load_maxrep 1 -load_cst 1  -s_path $< -t_path $(subst .s,.s,$<) >$@

clean : 
	rm -f $(TARGETS)
