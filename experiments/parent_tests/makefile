MS_EXEC   = ../../fast_ms/bin/matching_stats_parallel.x
SCRIPT = mstat_fast_ms.py
IDIR   = ../datasets/big_paper3

INPUT         = $(notdir $(subst .s,,$(wildcard $(IDIR)/rep*.s) $(wildcard $(IDIR)/rnd*.s)))
STATS_TARGETS = $(foreach opts,stats,$(addsuffix .$(opts).csv,$(INPUT)))
TIME_TARGETS  = $(foreach opts,pseq lca,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS          = --load_cst --repeat 3 --time_usage

rank_timing.csv : $(TIME_TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.lca.csv : $(IDIR)/%.s $(MS_EXEC)
	$(SCRIPT) --label lca  --lca_parents $(OPTS) --output $@ $< $(subst .s,,$<).t

%.pseq.csv : $(IDIR)/%.s $(MS_EXEC)
	$(SCRIPT) --label pseq               $(OPTS) --output $@ $< $(subst .s,,$<).t


clean : 
	rm -f  rank_timing.csv $(TARGETS)
