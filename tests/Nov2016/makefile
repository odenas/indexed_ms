
FAST_MS          = ../fd_ms
SLOW_MS          = ../slow_ms.py
BDIR             = input_data
LENS             = 400 800 1600 3200 6400 12800 25600
REP              = 15
PROTEIN_PREFIXES = $(addprefix proteins.50MB,$(foreach LS,$(LENS),$(foreach LT,$(LENS),$(LT)_$(LS))))
TARGETS          = $(addsuffix _protein_benchmark.csv,nonlazy lazy ohleb sada)

.PHONY: test clean random_input

all : $(TARGETS)


ohleb_protein_benchmark.csv :  ../script_benchmark.py $(FAST_MS) 
	python $< --repeat $(REP) --output $@ $(BDIR) --fast_prg $(FAST_MS) $(PROTEIN_PREFIXES)

sada_protein_benchmark.csv : ../script_benchmark.py $(FAST_MS)
	python $< --repeat $(REP) --output $@ $(BDIR) --fast_prg $(FAST_MS) --sada $(PROTEIN_PREFIXES)

lazy_protein_benchmark.csv :  ../script_benchmark.py $(FAST_MS)
	python $< --repeat $(REP) --lazy_wl --output $@ $(BDIR) --fast_prg $(FAST_MS) $(PROTEIN_PREFIXES)

nonlazy_protein_benchmark.csv : ../script_benchmark.py $(FAST_MS)
	python $< --repeat $(REP) --output $@ $(BDIR) --fast_prg $(FAST_MS) $(PROTEIN_PREFIXES)

clean  :
	rm -f $(TARGETS)


