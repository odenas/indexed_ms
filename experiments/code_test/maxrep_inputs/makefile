
DUMP_CST        = ../../../linux_build/dump_cst/main
DUMP_MAXREP     = ../../../linux_build/dump_maxrep/main

GENERATE_SCRIPT = ../../datasets/generate_input3.py
OPTS = --rep_mut 10 --rep_blocks 100 --sim_mut 200 --seed 42
SS   = 20

ALL_INPUT_PREF = $(foreach ts,10 11 12,rnd_$(SS)s_dis_$(ts)t_abcd)
ELST_INPUT_PREF = rnd_$(SS)s_dis_$(TS)t_$(ALPS)


INPUT_STR  = $(addsuffix .s,$(ALL_INPUT_PREF)) $(addsuffix .t,$(ALL_INPUT_PREF))
INPUT_CST  = $(addsuffix .s.fwd.stree,$(ALL_INPUT_PREF))
INPUT_MREP = $(addsuffix .s.rev.maxrep,$(ALL_INPUT_PREF))
INPUT_ELST = $(addsuffix .s.elst.00.bin,$(ALL_INPUT_PREF))
ALL_INPUT  = $(INPUT_STR) $(INPUT_CST) $(INPUT_MREP) $(INPUT_ELST)

.SECONDARY : $(ALL_INPUT)


#all : $(INPUT_ELST)
all : rnd_20s_dis_10t_abcd.s.elst.00.bin

rnd_%_abcd.s : $(GENERATE_SCRIPT)
	python $< rnd_dis abcd $(SS) $(subst t_abcd.s,,$(subst rnd_$(SS)s_dis_,,$@)) $(OPTS)

%.fwd.stree : % $(DUMP_CST)
	$(DUMP_CST) -s_path $< 

%.rev.maxrep : % %.fwd.stree $(DUMP_MAXREP)
	$(DUMP_MAXREP) -load_cst 1 -s_path $<

%.elst.00.bin : % %.rev.maxrep  ../../../linux_build/dump_edge_lists/main 
	../../../linux_build/dump_edge_lists/main -check 1 -s_path $<  -out_path $(subst .00.bin,,$@)

clean : 
	rm -f $(ALL_INPUT) *bin *stree

