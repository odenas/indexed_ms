
BDIR        = input_data
TARGETS     = $(addsuffix .threads.csv,$(addprefix results.,1 2 4 8 16))


.PHONY: all clean


all: $(TARGETS)

results.%.threads.csv : script_benchmark.py
	python $<   --load_cst  --nthreads $(subst .threads.csv,,$(subst results.,,$@)) --output $@  --repeat 5  $(BDIR) rnd_100Ms_5Mt 


clean  :
	rm -f $(TARGETS)

