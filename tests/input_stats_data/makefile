include ../global.makefile
SCRIPT  = input_stats.py
INPUT   = $(notdir $(subst .s,,$(wildcard $(IDIR)/rep*s) $(wildcard $(IDIR)/rnd*s)))
TARGETS = $(addsuffix .csv,$(INPUT))


stats.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.csv : $(SCRIPT) $(IDIR)/%.s $(IDIR)/%.t
	python $^ --load_maxrep --load_cst >$@

clean :
	rm -f $(TARGETS) stats.csv
