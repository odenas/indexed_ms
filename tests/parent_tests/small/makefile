
SCRIPT = ../../fast_ms.py
INPUT  = $(notdir $(subst .s,,$(wildcard input_data/rep*s) $(wildcard input_data/rnd*s)))
TARGETS= $(foreach opts,pseq lca,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS   = --load_cst --repeat 15 --space_usage --time_usage

#SANDBOX_TARGETS = rep_10000000s_sim_2000000t_abcd_sim1000.sandbox.csv rep_10000000s_sim_2000000t_abcdefghijklmnopqrst_sim1000.sandbox.csv
SANDBOX_TARGETS = rep_10000000s_sim_2000000t_abcd_sim1000.sandbox.csv

all : rank_timing.csv $(SANDBOX_TARGETS)

rank_timing.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.pseq.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --label pseq  $(OPTS) --output $@ $< $(subst .s,,$<).t

%.lca.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --label lca  --lca_parents $(OPTS) --output $@ $< $(subst .s,,$<).t

%.sandbox.csv : ../sandbox_timing.csv
	ln -s $< .

clean : 
	rm -f  rank_timing.csv $(TARGETS)
