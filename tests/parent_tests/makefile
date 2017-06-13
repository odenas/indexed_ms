
SCRIPT = ../fast_ms.py
INPUT  = $(notdir $(subst .s,,$(wildcard input_data/rep*s) $(wildcard input_data/rnd*s)))
TARGETS= $(foreach opts,pseq lca,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS   = --load_cst --repeat 15 --space_usage --time_usage

#SANDBOX_TARGETS = rep_10000000s_sim_2000000t_abcd_sim1000.sandbox.csv rep_10000000s_sim_2000000t_abcdefghijklmnopqrst_sim1000.sandbox.csv
SANDBOX_TARGETS = $(addsuffix .sandbox.csv,$(INPUT))

all : rank_timing.csv sandbox_timing.csv

sandbox_timing.csv : $(SANDBOX_TARGETS)
	awk -F "," 'NR==1 || $$1!="ntrial" {print $$0","(NR==1?"filename" : FILENAME)}' $^ >$@

rank_timing.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.pseq.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --label pseq  $(OPTS) --output $@ $< $(subst .s,,$<).t

%.lca.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --label lca  --lca_parents $(OPTS) --output $@ $< $(subst .s,,$<).t

%.sandbox.csv : ../xcode_bin/sanbox_parent_calls input_data/%.s
	$< -s_path input_data/$(subst .sandbox.csv,,$@).s -load_cst 1 >$@

clean : 
	rm -f  rank_timing.csv $(TARGETS)
