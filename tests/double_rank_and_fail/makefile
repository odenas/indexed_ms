
FAST_MS               = ../fast_ms.py
BDIR                  = input_data
REP                   = 10
PREFIXES              = $(notdir $(subst .s,,$(wildcard $(BDIR)/rnd*.s)))
SANDBOX_TARGETS       = $(addsuffix .csv,$(addprefix $(BDIR)/sandbox_timming_,$(notdir $(DATA_TARGETS))))


all: 
	"compile for each of the double_rank and double_rank_and_fail"

res.csv : double_rank.csv double_rank_and_fail.csv
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

double_rank.csv :  ./script_benchmark.py $(FAST_MS) $(DATA_TARGETS)
	python $< --load_cst --repeat $(REP) --output $@ $(BDIR) --label $(subst .csv,,$@) $(PREFIXES) --v >$@

double_rank_and_fail.csv :  ./script_benchmark.py $(FAST_MS) $(DATA_TARGETS)
	python $< --load_cst --repeat $(REP) --output $@ $(BDIR) --label $(subst .csv,,$@) $(PREFIXES) --v >$@

