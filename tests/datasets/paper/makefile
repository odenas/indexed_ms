
DUMP_CST        = ../../xcode_bin/dump_cst
DUMP_MAXREP     = ../../xcode_bin/dump_maxrep
GENERATE_SCRIPT = ../generate_input.py
KZEROS = 000
MZEROS = $(KZEROS)$(KZEROS)
GZEROS = $(MZEROS)$(KZEROS)


SMALL_ALPHABET = abcd
LARGE_ALPHABET = abcdefghijklmnopqrst

LENS           = 10
LENT           = 1
SIZE_PATTERN   = $(LENS)Ms_$(LENT)Mt
RND_INPUT_PREF = $(addprefix rnd_$(SIZE_PATTERN)_,$(SMALL_ALPHABET) $(LARGE_ALPHABET))
REP_INPUT_PREF = $(addprefix rep_$(SIZE_PATTERN)_,$(SMALL_ALPHABET) $(LARGE_ALPHABET))
MUT_INPUT_PREF = $(addprefix mut_$(SIZE_PATTERN)_1000_,$(SMALL_ALPHABET) $(LARGE_ALPHABET))
ALL_INPUT_PREF = $(MUT_INPUT_PREF) $(RND_INPUT_PREF) $(REP_INPUT_PREF)

INPUT_STR  = $(addsuffix .s,$(ALL_INPUT_PREF)) $(addsuffix .t,$(ALL_INPUT_PREF))
INPUT_CST  = $(addsuffix .s.fwd.stree,$(ALL_INPUT_PREF))
INPUT_MREP = $(addsuffix .s.rev.maxrep,$(ALL_INPUT_PREF))
ALL_INPUT  = $(INPUT_STR) $(INPUT_CST) $(INPUT_MREP)

.SECONDARY : $(ALL_INPUT)


all : $(INPUT_MREP)

mut_$(SIZE_PATTERN)_1000_%.s : $(GENERATE_SCRIPT)
	python $< mut $(subst .s,,$(subst mut_$(SIZE_PATTERN)_1000_,,$@)) \
		$@ $(LENS)$(MZEROS) $(subst .s,.t,$@) $(LENT)$(MZEROS) --mutation_period 10000

rnd_$(SIZE_PATTERN)_%.s : $(GENERATE_SCRIPT)
	python $< rnd $(subst .s,,$(subst rnd_$(SIZE_PATTERN)_,,$@)) \
		$@ $(LENS)$(MZEROS) $(subst .s,.t,$@) $(LENT)$(MZEROS)

rep_$(SIZE_PATTERN)_%.s : $(GENERATE_SCRIPT)
	python $< rep $(subst .s,,$(subst rep_$(SIZE_PATTERN)_,,$@)) \
		$@ $(LENS)$(MZEROS) $(subst .s,.t,$@) $(LENT)$(MZEROS) --seed_pool 2 --seed_length 10000

%.fwd.stree : % $(DUMP_CST)
	$(DUMP_CST) -s_path $< 

%.rev.maxrep : % %.fwd.stree $(DUMP_MAXREP)
	$(DUMP_MAXREP) -load_cst 1 -s_path $<

clean : 
	rm -f $(ALL_INPUT)
