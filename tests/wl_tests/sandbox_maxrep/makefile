
EXEC =../../../linux_build/sandbox_test_maxrep/main
SDIR = ../../datasets/small_paper3
BDIR = ../../datasets/big_paper3
R    = R CMD BATCH -q --vanilla --no-timing

TARGETS = $(addsuffix .csv,rep_100000000s_sim_500000t_abcd.s rep_1000000s_sim_200000t_abcd.s)

.PHONY : all rds

all : $(TARGETS)

rds : $(addsuffix .rds,maximal0_close0 maximal0_close1 maximal1_close0 maximal1_close1)


rep_100000000s_sim_500000t_abcd.s.csv : $(BDIR)/rep_100000000s_sim_500000t_abcd.s $(EXEC)
	$(EXEC) -load_maxrep 1 -load_cst 1  -s_path $< -t_path $(subst .s,.s,$<) >$@

rep_1000000s_sim_200000t_abcd.s.csv : $(SDIR)/rep_1000000s_sim_200000t_abcd.s $(EXEC)
	$(EXEC) -load_maxrep 1 -load_cst 1  -s_path $< -t_path $(subst .s,.s,$<) >$@

clean : 
	rm -f $(TARGETS)
