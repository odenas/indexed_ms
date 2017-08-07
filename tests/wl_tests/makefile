
include ../global.makefile

INPUT  = $(notdir $(subst .s,,$(wildcard $(IDIR)/rep*s) $(wildcard $(IDIR)/rnd*s)))
INPUT  = $(notdir $(subst .s,,$(wildcard $(IDIR)/r*st_sim1000.s)))
TARGETS= $(foreach opts,l0_f0_m0 l0_f1_m0 l0_f1_m1 l1_f0_m0 l1_f1_m0,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS   = --load_cst --repeat $(REPEAT) --space_usage --time_usage

rank_timing.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.l0_f0_m0.csv : $(IDIR)/%.s $(SCRIPT) $(BIN)
	python $(SCRIPT) --label l0_f0_m0  $(OPTS) --output $@ $< $(subst .s,,$<).t

#%.l0_f0_m1.csv : NA

%.l0_f1_m0.csv : $(IDIR)/%.s $(SCRIPT) $(BIN)
	python $(SCRIPT) --label l0_f1_m0 --rank_fail $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l0_f1_m1.csv : $(IDIR)/%.s $(SCRIPT) $(BIN)
	python $(SCRIPT) --label l0_f1_m1 --rank_fail --load_maxrep --use_maxrep $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l1_f0_m0.csv : $(IDIR)/%.s $(SCRIPT) $(BIN)
	python $(SCRIPT) --label l1_f0_m0 --lazy_wl $(OPTS) --output $@ $< $(subst .s,,$<).t

#%.l1_f0_m1.csv : NA

%.l1_f1_m0.csv : $(IDIR)/%.s $(SCRIPT)  $(BIN)
	python $(SCRIPT) --label l1_f1_m0 --lazy_wl --rank_fail $(OPTS) --output $@ $< $(subst .s,,$<).t

#%.l1_f1_m1.csv : NA

clean : 
	rm -f  rank_timing.csv $(TARGETS)
