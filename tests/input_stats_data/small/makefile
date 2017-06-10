
SCRIPT  = ../input_stats.py
INPUT   = $(notdir $(subst .s,,$(wildcard input_data/rep*s) $(wildcard input_data/rnd*s)))
TARGETS = $(addsuffix .csv,$(INPUT))


stats.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.csv : $(SCRIPT) input_data/%.s input_data/%.t
	python $^ --load_maxrep --load_cst >$@

clean :
	rm -f $(TARGETS) stats.csv
