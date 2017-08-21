
EP          = ../../xcode_bin/fd_ms
IDIR        = .
SCRIPT      = ../fast_ms.py
REPEAT      = 10
GENINP      = ./generate_input.py
DUMP_CST    = ../../../xcode_bin/dump_cst
DUMP_MAXREP = ../../../xcode_bin/dump_maxrep
ALP         = abcdefghijklmnopqrst
LEN_S       = $(addsuffix 000,100 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10000)
#LEN_S       = 80000 160000
LEN_T       = 2000000
MP          = 100
OPTS        = --load_cst --repeat $(REPEAT) --space_usage --time_usage --v --vv


INPUT_PREF         = $(foreach lens,$(LEN_S),rep_$(lens)s_sim_$(LEN_T)t_$(ALP)_sim$(MP))
INPUT_TARGETS      = $(addsuffix .s.rev.stree,$(INPUT_PREF)) $(addsuffix .s.fwd.stree,$(INPUT_PREF)) $(addsuffix .s,$(INPUT_PREF)) $(addsuffix .t,$(INPUT_PREF)) $(addsuffix .s.rev.maxrep,$(INPUT_PREF)) 
EXPERIMENT_TARGETS = $(addsuffix .l0_f1_m0.csv,$(INPUT_PREF)) $(addsuffix .l0_f1_m1.csv,$(INPUT_PREF))

.SECONDARY : $(EXPERIMENT_TARGETS) $(INPUT_TARGETS)

#all : $(addsuffix .s.rev.maxrep,$(INPUT_PREF))

stats.csv : $(EXPERIMENT_TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.l0_f1_m0.csv : $(IDIR)/%.s.rev.maxrep $(SCRIPT)
	python $(SCRIPT) --exec_path $(EP) --label l0_f1_m0 --rank_fail $(OPTS) --output $@ $(subst .rev.maxrep,,$<) $(subst .s.rev.maxrep,,$<).t

%.l0_f1_m1.csv : $(IDIR)/%.s.rev.maxrep $(SCRIPT)
	python $(SCRIPT) --exec_path $(EP) --label l0_f1_m1 --rank_fail $(OPTS) --output $@ --use_maxrep --load_maxrep $(subst .rev.maxrep,,$<) $(subst .s.rev.maxrep,,$<).t

rep_%_sim_$(LEN_T)t_$(ALP)_sim$(MP).s: $(GENINP) 
	python $< rep_sim $(ALP) $(subst s_sim_$(LEN_T)t_$(ALP)_sim$(MP).s,,$(subst rep_,,$@)) $(LEN_T) --mutation_period $(MP)

%.fwd.stree : % $(DUMP_CST)
	$(DUMP_CST) -s_path $< 

%.rev.maxrep : % %.fwd.stree $(DUMP_MAXREP)
	$(DUMP_MAXREP) -load_cst 1 -s_path $<

clean : 
	rm -f $(EXPERIMENT_TARGETS)
#rm -f $(INPUT_TARGETS) $(EXPERIMENT_TARGETS)
