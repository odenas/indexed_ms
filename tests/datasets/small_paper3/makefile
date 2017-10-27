
DUMP_CST        = ../../../linux_build/dump_cst/main
DUMP_MAXREP     = ../../../linux_build/dump_maxrep/main
DUMP_ELST       = ../../../linux_build/dump_edge_lists/main

GENERATE_SCRIPT = ../generate_input3.py

SS   = 1000000
TS   = 200000
# 100 blocks of length 10K, 1 mutation / 1K bases
OPTS = --rep_mut 10 --rep_blocks 100 --sim_mut 200 --seed 42

ALPS = abcd
ALPL = abcdefghijklmnopqrst

ALL_INPUT_PREF = $(foreach ts,rnd rep,$(foreach tt,sim dis,$(foreach alp,$(ALPS) $(ALPL),$(ts)_$(SS)s_$(tt)_$(TS)t_$(alp))))
ELST_INPUT_PREF = rep_$(SS)s_sim_$(TS)t_$(ALPS) rep_$(SS)s_sim_$(TS)t_$(ALPL)


INPUT_STR  = $(addsuffix .s,$(ALL_INPUT_PREF)) $(addsuffix .t,$(ALL_INPUT_PREF))
INPUT_CST  = $(addsuffix .s.fwd.stree,$(ALL_INPUT_PREF))
INPUT_MREP = $(addsuffix .s.rev.maxrep,$(ALL_INPUT_PREF))
INPUT_ELST = $(addsuffix .s.rev.elst.i00.bin,$(ELST_INPUT_PREF))
ALL_INPUT  = $(INPUT_STR) $(INPUT_CST) $(INPUT_MREP)

.SECONDARY : $(ALL_INPUT)


all : $(INPUT_MREP) $(INPUT_ELST)

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

%.rev.maxrep : % %.fwd.stree $(DUMP_MAXREP)
	$(DUMP_MAXREP) -load_cst 1 -s_path $<

%.rev.elst.i00.bin : % %.rev.maxrep $(DUMP_ELST)
	$(DUMP_ELST) -s_path $<  -load_cst 1 -load_maxrep 1 -sample_freq 10 -check 0

clean : 
	rm -f $(ALL_INPUT) \
		$(addsuffix .s.rev.stree,$(ALL_INPUT_PREF)) \
		$(foreach suff,i00 i01 i10 i11 l00 l01,$(addsuffix .s.rev.elst.$(suff).bin,$(ELST_INPUT_PREF)))
