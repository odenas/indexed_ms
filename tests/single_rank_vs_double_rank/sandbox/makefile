
SCRIPT = ../../xcode_bin/single_vs_double_rank
INPUT  = $(notdir $(subst .s,,$(wildcard ../input_data/rep*s) $(wildcard ../input_data/rnd*s)))
TARGETS= $(addsuffix .csv,$(INPUT))


rank_timing.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.csv : ../input_data/%.s $(SCRIPT) 
	$(SCRIPT) -load_cst 1 -s_path $< >$@

clean : 
	rm -f  rank_timing.csv $(TARGETS)
