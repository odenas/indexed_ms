
DUMP_CST        = ../../../../fast_ms/bin/dump_ms_data_structures.x
GENERATE_SCRIPT = ../../../datasets/generate_input3.py

# SLEN 100M, TLEN 500K
SS   = 100000000
TS   = 500000
# 10K blocks of length 10K, 1 mutation / 1K bases. 500 mutations on S to obtain T
OPTS = --rep_mut 10 --rep_blocks 10000 --sim_mut 500 --seed $(shell basename `pwd ` | cut -d _  -f 2 -)
#OPTS = --rep_mut 10 --rep_blocks 100 --sim_mut 5 --seed 42

# DEBUG #
#SS = 100000
#TS = 5000
## 1000 blocks of length 100, 1 mutation / 100 bases. 5 mutations on S to obtain T
#OPTS = --rep_mut 1 --rep_blocks 100 --sim_mut 5

ALPS = abcd
ALPL = abcdefghijklmnopqrst

ALL_INPUT_PREF  = $(foreach ts,rnd rep,$(foreach tt,sim dis,$(foreach alp,$(ALPS) $(ALPL),$(ts)_$(SS)s_$(tt)_$(TS)t_$(alp))))
INPUT_STR       = $(addsuffix .s,$(ALL_INPUT_PREF)) $(addsuffix .t,$(ALL_INPUT_PREF))
INPUT_CST       = $(addsuffix .s.fwd.stree,$(ALL_INPUT_PREF))
ALL_INPUT       = $(INPUT_STR) $(INPUT_CST)

.SECONDARY : $(ALL_INPUT)


all : $(INPUT_CST)

rnd_%_$(ALPS).s : $(GENERATE_SCRIPT)
	python $< $(subst _$(TS)t_$(ALPS).s,,$(subst $(SS)s_,,$@)) $(ALPS) $(SS) $(TS) $(OPTS)

rnd_%_$(ALPL).s : $(GENERATE_SCRIPT)
	python $< $(subst _$(TS)t_$(ALPL).s,,$(subst $(SS)s_,,$@)) $(ALPL) $(SS) $(TS) $(OPTS)

rep_%_$(ALPS).s : $(GENERATE_SCRIPT)
	python $< $(subst _$(TS)t_$(ALPS).s,,$(subst $(SS)s_,,$@)) $(ALPS) $(SS) $(TS) $(OPTS)

rep_%_$(ALPL).s : $(GENERATE_SCRIPT)
	python $< $(subst _$(TS)t_$(ALPL).s,,$(subst $(SS)s_,,$@)) $(ALPL) $(SS) $(TS) $(OPTS)

%.fwd.stree : % $(DUMP_CST)
	$(DUMP_CST) -s_path $< 

clean : 
	rm -f $(ALL_INPUT) \
		$(addsuffix .s.rev.stree,$(ALL_INPUT_PREF)) $(addsuffix .s.fwd.stree,$(ALL_INPUT_PREF)) \
		$(addsuffix .s.rev.maxrep,$(ALL_INPUT_PREF)) $(addsuffix .s.fwd.maxrep,$(ALL_INPUT_PREF))
