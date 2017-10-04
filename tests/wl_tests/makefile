
EXEC   = ../../linux_build/fd_ms/fd_ms
SCRIPT = ../fast_ms.py
IDIR   = ../datasets/big_paper2

INPUT  = $(notdir $(subst .s,,$(wildcard $(IDIR)/rep*s) $(wildcard $(IDIR)/rnd*s)))
TARGETS= $(foreach opts,l0_f0_m0 l0_f1_m0 l0_f1_m1 l1_f0_m0 l1_f1_m0,$(addsuffix .$(opts).csv,$(INPUT)))
REPEAT = 2
OPTS   = --load_cst --repeat $(REPEAT) --space_usage --time_usage


rank_timing.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.l0_f0_m0.csv : $(IDIR)/%.s $(SCRIPT) $(EXEC)
	python $(SCRIPT) --label l0_f0_m0  $(OPTS) --output $@ $< $(subst .s,,$<).t

#%.l0_f0_m1.csv : NA

%.l0_f1_m0.csv : $(IDIR)/%.s $(SCRIPT) $(EXEC)
	python $(SCRIPT) --label l0_f1_m0 --rank_fail $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l0_f1_m1.csv : $(IDIR)/%.s $(SCRIPT) $(EXEC)
	python $(SCRIPT) --label l0_f1_m1 --rank_fail --load_maxrep --use_maxrep $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l1_f0_m0.csv : $(IDIR)/%.s $(SCRIPT) $(EXEC)
	python $(SCRIPT) --label l1_f0_m0 --lazy_wl $(OPTS) --output $@ $< $(subst .s,,$<).t

#%.l1_f0_m1.csv : NA

%.l1_f1_m0.csv : $(IDIR)/%.s $(SCRIPT)  $(EXEC)
	python $(SCRIPT) --label l1_f1_m0 --lazy_wl --rank_fail $(OPTS) --output $@ $< $(subst .s,,$<).t

#%.l1_f1_m1.csv : NA

clean : 
	rm -f  rank_timing.csv $(TARGETS)
