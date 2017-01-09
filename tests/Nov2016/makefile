
FAST_MS          = ../fd_ms
BDIR             = input_data
REP              = 1
PROTEIN_PREFIXES = $(subst t.txt,,$(notdir $(wildcard $(BDIR)/*t.txt)))
TARGETS          = $(addsuffix .csv,$(addprefix $(BDIR)_,lazy ohleb sada))

.PHONY: test clean random_input

all : $(TARGETS)


$(BDIR)_ohleb.csv :  ../script_benchmark.py $(FAST_MS) 
	python $< --repeat $(REP) --output $@ $(BDIR) --fast_prg $(FAST_MS) --label ohleb $(PROTEIN_PREFIXES)

$(BDIR)_sada.csv : ../script_benchmark.py $(FAST_MS)
	python $< --repeat $(REP) --output $@ $(BDIR) --fast_prg $(FAST_MS) --sada --label sada $(PROTEIN_PREFIXES)

$(BDIR)_lazy.csv :  ../script_benchmark.py $(FAST_MS)
	python $< --repeat $(REP) --lazy_wl --output $@ $(BDIR) --fast_prg $(FAST_MS) --label lazy $(PROTEIN_PREFIXES)


clean  :
	rm -f $(TARGETS)


