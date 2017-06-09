
SCRIPT = ../../fast_ms.py
INPUT  = $(notdir $(subst .s,,$(wildcard input_data/rep*s) $(wildcard input_data/rnd*s)))
TARGETS= $(foreach opts,1 2 4 8 16,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS   = --load_cst --repeat 5 --space_usage --time_usage


.PHONY: clean

stats.csv: $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.1.csv : input_data/%.s $(SCRIPT)
	python $(SCRIPT) $(OPTS) --nthreads 1 --label 1 --output $@  $< $(subst .s,,$<).t

%.2.csv : input_data/%.s $(SCRIPT)
	python $(SCRIPT) $(OPTS) --nthreads 2 --label 2 --output $@  $< $(subst .s,,$<).t

%.4.csv : input_data/%.s $(SCRIPT)
	python $(SCRIPT) $(OPTS) --nthreads 4 --label 4 --output $@  $< $(subst .s,,$<).t

%.8.csv : input_data/%.s $(SCRIPT)
	python $(SCRIPT) $(OPTS) --nthreads 8 --label 8 --output $@  $< $(subst .s,,$<).t

%.16.csv : input_data/%.s $(SCRIPT)
	python $(SCRIPT) $(OPTS) --nthreads 16 --label 16 --output $@  $< $(subst .s,,$<).t

clean  :
	rm -f stats.csv $(TARGETS)
