
EXEC =../../../linux_build/sandbox_selectatdist/main
DDIR = ../../datasets/big_paper2
#INPUT = $(DDIR)/rep_1000000s_dis_200000t_abcd_sim1000.s

INPUT  = $(notdir $(subst .s,,$(wildcard $(DDIR)/rep*t_sim1000.s) $(wildcard $(DDIR)/rnd*t_sim1000.s)))
TARGETS = $(addsuffix .csv,$(INPUT))

all : $(TARGETS)

%.csv : $(DDIR)/%.s $(EXEC)
	$(EXEC) -load_maxrep 1 -load_cst 1  -s_path $< >$@

clean : 
	rm -f $(TARGETS)
