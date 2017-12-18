MS_EXEC   = ../../fast_ms/bin/matching_stats.x
ST_EXEC   = ../../fast_ms/bin/parent_seq_stats.x
SCRIPT = ../fast_ms.py
IDIR   = ../datasets/big_paper3

INPUT         = $(notdir $(subst .s,,$(wildcard $(IDIR)/rep*s) $(wildcard $(IDIR)/rnd*s)))
STATS_TARGETS = $(foreach opts,stats,$(addsuffix .$(opts).csv,$(INPUT)))
TIME_TARGETS  = $(foreach opts,pseq lca,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS          = --load_cst --repeat 3 --time_usage

stats.csv : $(STATS_TARGETS)
	awk -F "," 'NR==1 || $$1!="method"{print $$0","FILENAME}' $^ >$@

rank_timing.csv : $(TIME_TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.lca.csv : $(IDIR)/%.s $(SCRIPT) $(MS_EXEC)
	python $(SCRIPT) --label lca  --lca_parents $(OPTS) --output $@ $< $(subst .s,,$<).t

%.pseq.csv : $(IDIR)/%.s $(SCRIPT) $(MS_EXEC)
	python $(SCRIPT) --label pseq               $(OPTS) --output $@ $< $(subst .s,,$<).t

%.stats.csv : $(IDIR)/%.s  $(ST_EXEC)
	$(ST_EXEC) -s_path $< -t_path $(subst .s,.t,$<) -load_cst 1 >$@

clean : 
	rm -f  rank_timing.csv $(TARGETS)
