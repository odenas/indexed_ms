
#SCRIPT = script_benchmark.py
SCRIPT = ../fast_ms.py
INPUT  = $(notdir $(subst .s,,$(wildcard input_data/rep*s) $(wildcard input_data/rnd*s)))
TARGETS= $(addsuffix .csv,$(INPUT)) $(addsuffix .fail.csv,$(INPUT))
OPTS   = --load_cst --repeat 5 --space_usage --time_usage

rank_timing.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.fail.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --label double_rank_fail --rank_fail $(OPTS) --output $@ input_data/$(subst .fail.csv,,$@).s input_data/$(subst .fail.csv,,$@).t

%.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --label double_rank   $(OPTS) --output $@ input_data/$(subst .csv,,$@).s input_data/$(subst .csv,,$@).t

clean : 
	rm -f  rank_timing.csv $(TARGETS)
