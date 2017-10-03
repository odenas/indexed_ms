
DUMP_CST        = ../../xcode_bin/dump_cst
DUMP_MAXREP     = ../../../linux_build/dump_maxrep/main
GENERATE_SCRIPT = ../generate_input2.py
MZEROS = 000000

SS   = 100000000
TS   = 500000
MP   = 1000

ALPS = abcd
ALPL = abcdefghijklmnopqrst

ALL_INPUT_PREF = $(foreach ts,rnd rep,$(foreach tt,sim dis,$(foreach alp,$(ALPS) $(ALPL),$(ts)_$(SS)s_$(tt)_$(TS)t_$(alp)_sim$(MP))))


INPUT_STR  = $(addsuffix .s,$(ALL_INPUT_PREF)) $(addsuffix .t,$(ALL_INPUT_PREF))
INPUT_CST  = $(addsuffix .s.fwd.stree,$(ALL_INPUT_PREF))
INPUT_MREP = $(addsuffix .s.rev.maxrep,$(ALL_INPUT_PREF))
ALL_INPUT  = $(INPUT_STR) $(INPUT_CST) $(INPUT_MREP)

.SECONDARY : $(ALL_INPUT)


all : $(INPUT_MREP)

rnd_%_$(ALPS)_sim$(MP).s : $(GENERATE_SCRIPT)
	python $< $(subst _$(TS)t_$(ALPS)_sim$(MP).s,,$(subst $(SS)s_,,$@)) $(ALPS) $(SS) $(TS) --mutation_period $(MP) --seed_pool 2 --seed_length 10000

rnd_%_$(ALPL)_sim$(MP).s : $(GENERATE_SCRIPT)
	python $< $(subst _$(TS)t_$(ALPL)_sim$(MP).s,,$(subst $(SS)s_,,$@)) $(ALPL) $(SS) $(TS) --mutation_period $(MP) --seed_pool 2 --seed_length 10000

rep_%_$(ALPS)_sim$(MP).s : $(GENERATE_SCRIPT)
	python $< $(subst _$(TS)t_$(ALPS)_sim$(MP).s,,$(subst $(SS)s_,,$@)) $(ALPS) $(SS) $(TS) --mutation_period $(MP) --seed_pool 2 --seed_length 10000

rep_%_$(ALPL)_sim$(MP).s : $(GENERATE_SCRIPT)
	python $< $(subst _$(TS)t_$(ALPL)_sim$(MP).s,,$(subst $(SS)s_,,$@)) $(ALPL) $(SS) $(TS) --mutation_period $(MP) --seed_pool 2 --seed_length 10000

%.fwd.stree : % $(DUMP_CST)
	$(DUMP_CST) -s_path $< 

%.rev.maxrep : % %.fwd.stree $(DUMP_MAXREP)
	$(DUMP_MAXREP) -load_cst 1 -s_path $<

clean : 
	rm -f $(ALL_INPUT)
