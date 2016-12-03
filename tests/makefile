
FAST_MS = ./fd_ms
SLOW_MS = slow_ms.py
BDIR    = input_data
PROTEIN_DSETS = $(addprefix proteins.50MB1000_,1000 5000 10000 20000 40000 80000 160000 320000)
N = 10

.PHONY: test clean input

protein_benchmark.csv :  script_benchmark.py $(FAST_MS)
	python $< --output $@ $(BDIR) $(PROTEIN_DSETS)

test: script_tests.py $(SLOW_MS) $(FAST_MS)  | $(BDIR)
	python $< --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) $| $(shell python -c "print ' '.join(map(str, range($(N))))")

random_input : generate_input.py | $(BDIR)
	python $< random --l 400  --L 800 --base_dir $| --n $(N)

clean  :
	rm -f a.txt

