
FAST_MS          = ./fd_ms
SLOW_MS          = slow_ms.py
BDIR             = input_data
RANDOM_ALPHABET  = ab
PROTEIN_PREFIXES = $(addprefix proteins.50MB1000_,1000 5000 10000 20000 40000 80000 160000 320000)
RANDOM_S_LENS    = 2 4 8 16 32 64 128 256 512 1024 2048
RANDOM_T_LEN     = 200
RANDOM_PREFIXES  = $(addprefix $(RANDOM_ALPHABET)$(RANDOM_T_LEN)_,$(RANDOM_S_LENS))

.PHONY: test clean random_input

protein_benchmark.csv :  script_benchmark.py $(FAST_MS)
	python $< --output $@ $(BDIR) $(PROTEIN_PREFIXES)

test: script_tests.py $(SLOW_MS) $(FAST_MS)  | $(BDIR)
	python $< --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) $| $(RANDOM_PREFIXES)

random_input : generate_input.py | $(BDIR)
	for sl in $(RANDOM_S_LENS); do python $< random $(RANDOM_ALPHABET) --len_t $(RANDOM_T_LEN)  --len_s $$sl --base_dir $|; done

clean  :
	rm -f a.txt

