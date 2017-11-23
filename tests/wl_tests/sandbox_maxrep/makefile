
EXEC = ../../../fast_ms/bin/sandbox_maxrep.x
SDIR = ../../datasets/small_paper3
BDIR = ../../datasets/big_paper3
R    = R CMD BATCH -q --vanilla --no-timing

TARGETS = $(addsuffix .csv, rep_1000000s_sim_200000t_abcd.s rep_1000000s_sim_200000t_abcdefghijklmnopqrst.s) \
		  $(addsuffix .csv, rep_100000000s_sim_500000t_abcd.s rep_100000000s_sim_500000t_abcdefghijklmnopqrst.s)


.PHONY : all rds

all : $(TARGETS)

rds : $(addsuffix .rds,maximal0_close0 maximal0_close1 maximal1_close0 maximal1_close1)


rep_100000000s_sim_500000t_abcd.s.csv : $(BDIR)/rep_100000000s_sim_500000t_abcd.s $(EXEC)
	$(EXEC) -load_maxrep 1 -load_cst 1 -repeat 5 -s_path $< >$@

rep_100000000s_sim_500000t_abcdefghijklmnopqrst.s.csv : $(BDIR)/rep_100000000s_sim_500000t_abcdefghijklmnopqrst.s $(EXEC)
	$(EXEC) -load_maxrep 1 -load_cst 1 -repeat 2  -s_path $< >$@

rep_1000000s_sim_200000t_abcd.s.csv : $(SDIR)/rep_1000000s_sim_200000t_abcd.s $(EXEC)
	$(EXEC) -load_maxrep 1 -load_cst 1 -repeat 5 -s_path $< >$@

rep_1000000s_sim_200000t_abcdefghijklmnopqrst.s.csv : $(SDIR)/rep_1000000s_sim_200000t_abcdefghijklmnopqrst.s $(EXEC)
	$(EXEC) -load_maxrep 1 -load_cst 1 -repeat 2 -s_path $< >$@

clean : 
	rm -f $(TARGETS)
