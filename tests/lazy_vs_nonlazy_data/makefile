
FAST_MS  = ../fast_ms.py
BDIR     = input_data
REP      = 1
PREFIXES = $(subst t.txt,,$(notdir $(wildcard $(BDIR)/*t.txt)))
TARGETS  = $(addsuffix .csv,$(addprefix $(BDIR)_,lazy nonlazy))

## dataset
TLEN = 20000
SLEN = 30000 40000 50000
RANDOM_DATA_TARGETS = $(addsuffix _$(TLEN).s,$(addprefix random_,$(SLEN)))
MAX_SLEN = 50000
MUT_PERIODS = 10 80 320 1280
MUTATION_DATA_TARGETS = $(addsuffix .s,$(addprefix mutation_$(MAX_SLEN)_$(TLEN)_,$(MUT_PERIODS)))


.PHONY: all clean

all : $(RANDOM_DATA_TARGETS) $(MUTATION_DATA_TARGETS)


#$(BDIR)_nonlazy.csv :  ./script_benchmark.py $(FAST_MS) 
#	python $< --repeat $(REP) --output $@ $(BDIR) --fast_prg $(FAST_MS) --label non_lazy --ms_progress 10 $(PREFIXES)

#$(BDIR)_lazy.csv :  ./script_benchmark.py $(FAST_MS)
#	python $< --repeat $(REP) --output $@ $(BDIR) --fast_prg $(FAST_MS) --label lazy --lazy_wl --ms_progress 10 $(PREFIXES)


mutation_$(MAX_SLEN)_$(TLEN)_%.s : 
	python ../generate_input.py random abcdefghijklmnopqrst \
		$@ $(MAX_SLEN) $(subst .s,.t,$@) $(TLEN) \
		--mutation_period $(subst .s,,$(subst mutation_$(MAX_SLEN)_$(TLEN)_,,$@))

random_%_$(TLEN).s : 
	python ../generate_input.py random abcdefghijklmnopqrst \
		$@ $(subst _$(TLEN).s,,$(subst random_,,$@)) $(subst .s,.t,$@) $(TLEN)


clean  :
	rm -f $(TARGETS)


