
SCRIPT = script_benchmark.py
INPUT  = $(notdir $(subst .s,,$(wildcard input_data/rep*s) $(wildcard input_data/rnd*s)))
TARGETS= $(addsuffix .csv,$(INPUT)) $(addsuffix .fail.csv,$(INPUT))
OPTS   = --load_cst --repeat 5 --space_usage --time_usage


res.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

double_rank_fail.csv :  ./script_benchmark.py $(addsuffix .s,$(addprefix input_data/,$(INPUT)))
	python $< --load_cst --repeat $(REPEAT) --output $@ --label $(subst .csv,,$@) input_data $(INPUT)

double_rank_fail_maxrep.csv :  ./script_benchmark.py $(addsuffix .s,$(addprefix input_data/,$(INPUT)))
	python $< --load_cst --use_maxrep --load_maxrep --repeat $(REPEAT) --output $@ --label $(subst .csv,,$@) input_data $(INPUT)

%.maxrep.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --label double_rank_maxrep $(OPTS) --output $@ input_data/$(subst .csv,,$@).s input_data/$(subst .csv,,$@).t

%.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --label double_rank $(OPTS) --output $@ input_data/$(subst .csv,,$@).s input_data/$(subst .csv,,$@).t

clean  : 
	rm -f res.csv double_rank_fail.csv double_rank_fail_maxrep.csv
