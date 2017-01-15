
include ../makefile.paths

BDIR     = input_data
REP      = 1
PREFIXES = $(subst t.txt,,$(notdir $(wildcard $(BDIR)/*t.txt)))
TARGETS  = $(addsuffix .csv,$(addprefix $(BDIR)_,lazy nonlazy))

.PHONY: all clean

all : $(TARGETS)


$(BDIR)_nonlazy.csv :  ../script_benchmark.py $(FAST_MS) 
	python $< --repeat $(REP) --output $@ $(BDIR) --fast_prg $(FAST_MS) --label non_lazy --ms_progress 10 $(PREFIXES)

$(BDIR)_lazy.csv :  ../script_benchmark.py $(FAST_MS)
	python $< --repeat $(REP) --output $@ $(BDIR) --fast_prg $(FAST_MS) --label lazy --lazy_wl --ms_progress 10 $(PREFIXES)


clean  :
	rm -f $(TARGETS)


