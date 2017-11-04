
DDIR   = ../datasets/big_paper2
SCRIPT = ../../linux_build/sandbox_parent_calls/main

INPUT  = $(notdir $(subst .s,,$(wildcard $(DDIR)/rep*s) $(wildcard $(DDIR)/rnd*s)))
TARGETS= $(foreach opts,pseq lca,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS   = --load_cst --repeat 5 --space_usage --time_usage


all : rank_timing.csv

#sandbox_timing.csv : $(SANDBOX_TARGETS)
#	awk -F "," 'NR==1 || $$1!="ntrial" {print $$0","(NR==1?"filename" : FILENAME)}' $^ >$@

rank_timing.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.pseq.csv : $(DDIR)/%.s $(SCRIPT) 
	python $(SCRIPT) --label pseq  $(OPTS) --output $@ $< $(subst .s,,$<).t

%.lca.csv : $(DDIR)/%.s $(SCRIPT) 
	python $(SCRIPT) --label lca  --lca_parents $(OPTS) --output $@ $< $(subst .s,,$<).t

clean : 
	rm -f  rank_timing.csv $(TARGETS)
