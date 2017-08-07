
DUMP_CST        = ../../xcode_bin/dump_cst
DUMP_MAXREP     = ../../xcode_bin/dump_maxrep
GENERATE_SCRIPT = ../generate_input.py
KZEROS = 000
MZEROS = $(KZEROS)$(KZEROS)
GZEROS = $(MZEROS)$(KZEROS)
ALP = abcdefghijklmnopqrst
MUT_PERIODS     = 10 20 30 40 50 60 70 80 90 100 1000 10000 100000 500000
MUT_INPUT_PREF  = $(addprefix mut_100Ms_5Mt_,$(MUT_PERIODS)) $(addprefix mut_1Gs_5Mt_,$(MUT_PERIODS))
RND_INPUT_PREF  = $(addprefix rnd_,1Ms_5Mt 100Ms_500Kt 100Ms_50Kt 100Ms_5Mt 1Gs_5Mt 1Gs_500Kt)
REP_INPUT_PREF  = $(addprefix rep_,100Ms_500Kt 100Ms_5Mt)
ALL_INPUT_PREF  = $(MUT_INPUT_PREF) $(RND_INPUT_PREF) $(REP_INPUT_PREF)

INPUT_STR  = $(addsuffix .s,$(ALL_INPUT_PREF))
INPUT_CST  = $(addsuffix .s.fwd.stree,$(ALL_INPUT_PREF))
INPUT_MREP = $(addsuffix .s.rev.maxrep,$(ALL_INPUT_PREF))
ALL_INPUT  = $(INPUT_STR) $(INPUT_CST) $(INPUT_MREP)

.SECONDARY : $(ALL_INPUT)


all : $(INPUT_MREP)

mut_100Ms_5Mt_%.s : $(GENERATE_SCRIPT)
	python $< mut $(ALP) $@ 100$(MZEROS) $(subst .s,.t,$@) 5$(MZEROS) --mutation_period $(subst .s,,$(subst mut_100Ms_5Mt_,,$@))

mut_1Gs_5Mt_%.s : $(GENERATE_SCRIPT)
	python $< mut $(ALP) $@ 1$(GZEROS) $(subst .s,.t,$@) 5$(MZEROS) --mutation_period $(subst .s,,$(subst mut_1Gs_5Mt_,,$@))

rnd_100Ms_5Mt_%.s : $(GENERATE_SCRIPT)
	python $< rnd $(ALP) $@ 100$(MZEROS) $(subst .s,.t,$@) 5$(MZEROS) --mutation_period $(subst .s,,$(subst mut_100Ms_5Mt_,,$@))

rnd_1Gs_5Mt.s : $(GENERATE_SCRIPT)
	python $< rnd $(ALP) $@ 1$(GZEROS) $(subst .s,.t,$@) 5$(MZEROS)

rnd_100Ms_5Mt.s : $(GENERATE_SCRIPT)
	python $< rnd $(ALP) $@ 100$(MZEROS) $(subst .s,.t,$@) 5$(MZEROS)

rnd_1Ms_5Mt.s : $(GENERATE_SCRIPT)
	python $< rnd $(ALP) $@ 1$(MZEROS) $(subst .s,.t,$@) 5$(MZEROS)

rnd_1Gs_500Kt.s : $(GENERATE_SCRIPT)
	python $< rnd $(ALP) $@ 1$(GZEROS) $(subst .s,.t,$@) 500$(KZEROS)

rep_100Ms_500Kt.s : $(GENERATE_SCRIPT)
	python $< rep $(ALP) $@ 100$(MZEROS) $(subst .s,.t,$@) 500$(KZEROS)

rep_100Ms_5Mt.s : $(GENERATE_SCRIPT)
	python $< rep $(ALP) $@ 100$(MZEROS) $(subst .s,.t,$@) 5$(MZEROS)

%.fwd.stree : % $(DUMP_CST)
	$(DUMP_CST) $< 

%.rev.maxrep : % %.fwd.stree $(DUMP_MAXREP)
	$(DUMP_MAXREP) -load_cst 1 -s_path $<
