
FAST_MS               = ../fast_ms.py
BDIR                  = input_data
REP                   = 8
#PREFIXES              = $(notdir $(subst .s,,$(wildcard $(BDIR)/rnd_*.s)))
PREFIXES              = $(notdir $(subst .s,,$(wildcard $(BDIR)/rnd_1M*.s)))
SANDBOX_TARGETS       = $(addsuffix .csv,$(addprefix $(BDIR)/sandbox_timming_,$(notdir $(DATA_TARGETS))))


res.csv : double_rank.csv double_rank_maxrep.csv
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

double_rank.csv :  ./script_benchmark.py $(FAST_MS) $(DATA_TARGETS)
	python $< --load_cst --repeat $(REP) --output $@ $(BDIR) --label $(subst .csv,,$@) $(PREFIXES) --v

double_rank_maxrep.csv :  ./script_benchmark.py $(FAST_MS) $(DATA_TARGETS)
	python $< --use_maxrep --load_cst --repeat $(REP) --output $@ $(BDIR) --label $(subst .csv,,$@) $(PREFIXES) --v

clean  : 
	rm -f res.csv double_rank_maxrep.csv double_rank.csv
