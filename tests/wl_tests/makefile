
EXEC   = ../../fast_ms/bin/matching_stats.opt.x
SCRIPT = ../fast_ms.py
IDIR   = ../datasets/big_paper3

INPUT  = $(notdir $(subst .s,,$(wildcard $(IDIR)/rep*s) $(wildcard $(IDIR)/rnd*s)))
LABELS = l0_f1_r1_m1 l0_f0_r1_m1 l0_f1_r1_m2 l0_f0_r1_m2 \
		 l1_f1_r1_m0 l1_f0_r1_m0 l1_f0_r0_m0 l0_f1_r1_m0 l0_f0_r1_m0 l0_f0_r0_m0
TARGETS= $(foreach lab,$(LABELS),$(addsuffix .$(lab).csv,$(INPUT)))
REPEAT = 5
OPTS   = --load_cst --load_maxrep --repeat $(REPEAT) --time_usage --lca_parents 


rank_timing.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

# both types of maxrep strategies: doublerank, fail, no lazy. m1 is vanilla, m2 is rank and check
# l0_f1_r1_m1
# l0_f0_r1_m1
# l0_f1_r1_m2
# l0_f0_r1_m2
%.l0_f1_r1_m1.csv : $(IDIR)/%.s $(SCRIPT)  $(EXEC)
	python $(SCRIPT) --label l0_f1_r1_m1 --rank_fail --double_rank --use_maxrep_vanilla $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l0_f0_r1_m1.csv : $(IDIR)/%.s $(SCRIPT)  $(EXEC)
	python $(SCRIPT) --label l0_f0_r1_m1             --double_rank --use_maxrep_vanilla $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l0_f1_r1_m2.csv : $(IDIR)/%.s $(SCRIPT)  $(EXEC)
	python $(SCRIPT) --label l0_f1_r1_m1 --rank_fail --double_rank --use_maxrep_rc      $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l0_f0_r1_m2.csv : $(IDIR)/%.s $(SCRIPT)  $(EXEC)
	python $(SCRIPT) --label l0_f0_r1_m1             --double_rank --use_maxrep_rc      $(OPTS) --output $@ $< $(subst .s,,$<).t


# single rank vs. double rank exclude maximal info: lazy_fail_doublerank
# l1_f1_r1_m0
# l1_f1_r0_m0 x
# l1_f0_r1_m0
# l1_f0_r0_m0
# l0_f1_r1_m0 
# l0_f1_r0_m0 x
# l0_f0_r1_m0
# l0_f0_r0_m0
%.l1_f1_r1_m0.csv : $(IDIR)/%.s $(SCRIPT)  $(EXEC)
	python $(SCRIPT) --label l1_f1_r1_m0 --lazy_wl --rank_fail --double_rank $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l1_f0_r1_m0.csv : $(IDIR)/%.s $(SCRIPT)  $(EXEC)
	python $(SCRIPT) --label l1_f0_r1_m0 --lazy                --double_rank $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l1_f0_r0_m0.csv : $(IDIR)/%.s $(SCRIPT)  $(EXEC)
	python $(SCRIPT) --label l1_f0_r0_m0 --lazy_wl                           $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l0_f1_r1_m0.csv : $(IDIR)/%.s $(SCRIPT)  $(EXEC)
	python $(SCRIPT) --label l0_f1_r1_m0           --rank_fail --double_rank $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l0_f0_r1_m0.csv : $(IDIR)/%.s $(SCRIPT)  $(EXEC)
	python $(SCRIPT) --label l0_f0_r1_m0                       --double_rank $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l0_f0_r0_m0.csv : $(IDIR)/%.s $(SCRIPT)  $(EXEC)
	python $(SCRIPT) --label l0_f0_r0_m0                                     $(OPTS) --output $@ $< $(subst .s,,$<).t


clean : 
	rm -f  rank_timing.csv $(TARGETS)
