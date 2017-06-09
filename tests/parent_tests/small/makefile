
SCRIPT = ../../fast_ms.py
INPUT  = $(notdir $(subst .s,,$(wildcard input_data/rep*s) $(wildcard input_data/rnd*s)))
TARGETS= $(foreach opts,pseq lca,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS   = --load_cst --repeat 5 --space_usage --time_usage

rank_timing.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.pseq.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --label pseq  $(OPTS) --output $@ $< $(subst .s,,$<).t

%.lca.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --label lca  --lca_parents $(OPTS) --output $@ $< $(subst .s,,$<).t

clean : 
	rm -f  rank_timing.csv $(TARGETS)
