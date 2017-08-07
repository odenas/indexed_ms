
SCRIPT = ./fast_ms.py
IDIR   = ./
REPEAT = 3
EP     = ../../xcode_bin/fd_ms

INPUT  = $(notdir $(subst .s,,$(wildcard $(IDIR)/rep*s) $(wildcard $(IDIR)/rnd*s)))
INPUT  = $(notdir $(subst .s,,$(wildcard $(IDIR)/re*st_sim1000.s)))
INPUT  = rep_50000000s_sim_5000000t_abcdefghijklmnopqrst_sim1000
TARGETS= $(foreach opts,l0_f1_m0 l0_f1_m1,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS   = --load_cst --repeat $(REPEAT) --space_usage --time_usage --v --vv

#rank_and_check__uv single_rank__uv rank_and_check__has_wl single_rank__has_wl

rank_timing.$(notdir $(EP)).csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

#%.l0_f0_m0.csv : NA

#%.l0_f0_m1.csv : NA

%.l0_f1_m0.csv : $(IDIR)/%.s $(SCRIPT) $(EP)
	python $(SCRIPT) --exec_path $(EP) --label l0_f1_m0 --rank_fail $(OPTS) --output $@ $< $(subst .s,,$<).t

%.l0_f1_m1.csv : $(IDIR)/%.s $(SCRIPT) $(EP)
	python $(SCRIPT) --exec_path $(EP) --label l0_f1_m1 --rank_fail $(OPTS) --output $@ --use_maxrep --load_maxrep $< $(subst .s,,$<).t

#%.l1_f0_m0.csv : NA

#%.l1_f0_m1.csv : NA

#%.l1_f1_m0.csv : NA

#%.l1_f1_m1.csv : NA

clean : 
	rm -f $(TARGETS)
