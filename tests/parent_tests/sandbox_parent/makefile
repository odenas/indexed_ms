
DDIR = ../../datasets/small_paper3
EXEC = ../../../linux_build/sandbox_parent_calls/main
R    = R CMD BATCH -q --vanilla --no-timing
#INPUT = $(DDIR)/rep_1000000s_dis_200000t_abcd_sim1000.s
#                rep_100000000s_dis_500000t_abcdefghijklmnopqrst_sim1000

INPUT  = rep_1000000s_sim_200000t_abcdefghijklmnopqrst rep_1000000s_sim_200000t_abcd
TARGETS = $(addsuffix .csv,$(INPUT))


all.csv : $(TARGETS)
	awk -F "," 'NR==1 || $$1!="parent" {print $$0","(NR==1?"filename" : FILENAME)}' $^ >$@

%.csv : $(DDIR)/%.s $(EXEC)
	$(EXEC) -load_cst 1 -repeat 1  -s_path $< >$@

clean : 
	rm -f $(TARGETS)
