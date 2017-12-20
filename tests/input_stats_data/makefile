
EXEC    = ../../fast_ms/bin/input_stats.x
IDIR    = ../datasets/big_paper3
INPUT   = $(notdir $(subst .s,,$(wildcard $(IDIR)/rep*s) $(wildcard $(IDIR)/rnd*s)))
TARGETS = $(addsuffix .csv,$(INPUT))


stats.csv : $(TARGETS)
	awk -F "," 'NR==1 || $$1!="len_s"{print $$0","FILENAME}' $^ >$@

%.csv : $(IDIR)/%.s $(IDIR)/%.t $(EXEC) 
	$(EXEC) -s_path $<  -t_path $(subst .s,.t,$<) -load_maxrep  1 -load_cst 1 >$@

clean :
	rm -f $(TARGETS) stats.csv
