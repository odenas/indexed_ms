
LAZY_TIMMING_SCRIPT   = /Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/lazy_wl_performance
FAST_MS               = ../fast_ms.py
BDIR                  = input_data
REP                   = 1
## dataset
#TLEN                  = 200000
#SLEN                  = $(addsuffix 000000,32 64 128 256 512)
#RANDOM_DATA_TARGETS   = $(addsuffix _$(TLEN).s,$(addprefix $(BDIR)/random_,$(SLEN)))
#MAX_SLEN              = $(shell python -c "print max(map(int, '$(SLEN)'.split()))")
#MUT_PERIODS           = 10 100 1000 10000 100000
#MUTATION_DATA_TARGETS = $(addsuffix .s,$(addprefix $(BDIR)/mutation_$(MAX_SLEN)_$(TLEN)_,$(MUT_PERIODS)))
DATA_TARGETS          = $(wildcard $(BDIR)/*.s)
RND_TARGETS           = $(wildcard $(BDIR)/rnd*.s)
PREFIXES              = $(notdir $(subst .s,,$(DATA_TARGETS)))
SANDBOX_TARGETS       = $(addsuffix .csv,$(addprefix $(BDIR)/sandbox_timming_,$(notdir $(DATA_TARGETS))))


.PHONY: all clean


all: lazy_vs_nonlazy.csv sandbox_timing.csv

sandbox_timing.csv : $(RND_TARGETS)
	$(LAZY_TIMMING_SCRIPT) -s_path $< >$@

lazy_vs_nonlazy.csv : lazy.csv nonlazy.csv
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

lazy.csv :  ./script_benchmark.py $(FAST_MS) $(DATA_TARGETS)
	python $< --repeat $(REP) --output $@ $(BDIR) --label lazy --lazy_wl --runs_progress 5 --ms_progress 5 $(PREFIXES) --v

nonlazy.csv :  ./script_benchmark.py $(FAST_MS) $(DATA_TARGETS)
	python $< --repeat $(REP) --output $@ $(BDIR) --label nonlazy --runs_progress 5 --ms_progress 5 $(PREFIXES) --v


$(BDIR)/mutation_$(MAX_SLEN)_$(TLEN)_%.s : 
	python ../generate_input.py random abcdefghijklmnopqrst \
		$@ $(MAX_SLEN) $(subst .s,.t,$@) $(TLEN) \
		--mutation_period $(subst .s,,$(subst $(BDIR)/mutation_$(MAX_SLEN)_$(TLEN)_,,$@))

$(BDIR)/random_%_$(TLEN).s : 
	python ../generate_input.py random abcdefghijklmnopqrst \
		$@ $(subst _$(TLEN).s,,$(subst $(BDIR)/random_,,$@)) $(subst .s,.t,$@) $(TLEN)


clean  :
	rm -f $(subst .s,.t,$(DATA_TARGETS)) $(DATA_TARGETS)


