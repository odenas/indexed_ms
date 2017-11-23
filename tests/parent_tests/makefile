EXEC   = ../../linux_build/fd_ms/fd_ms
SCRIPT = ../fast_ms.py
IDIR   = ../datasets/big_paper3

INPUT  = $(notdir $(subst .s,,$(wildcard $(IDIR)/rep*s) $(wildcard $(IDIR)/rnd*s)))
TARGETS= $(foreach opts,lca pseq,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS   = --load_cst --repeat 5 --time_usage


rank_timing.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.lca.csv : $(IDIR)/%.s $(SCRIPT) $(EXEC)
	python $(SCRIPT) --label lca  $(OPTS) --output $@ $< $(subst .s,,$<).t


%.pseq.csv : $(IDIR)/%.s $(SCRIPT) $(EXEC)
	python $(SCRIPT) --label pseq --rank_fail $(OPTS) --output $@ $< $(subst .s,,$<).t

clean : 
	rm -f  rank_timing.csv $(TARGETS)
