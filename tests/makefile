
FAST_MS          = ./fd_ms
SLOW_MS          = slow_ms.py
BDIR             = debug_input_data
#PROTEIN_PREFIXES = $(addprefix proteins.50MB1000_,1000 5000 10000 20000 40000 80000 160000 320000)
PROTEIN_PREFIXES = $(addprefix proteins.50MB, 12000_1024000)

RANDOM_ALPHABET  = ab
RANDOM_BDIR      = input_data
RANDOM_S_LENS    = 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768
RANDOM_T_LEN     = 200
RANDOM_PREFIXES  = $(addprefix $(RANDOM_ALPHABET)$(RANDOM_T_LEN)_,$(RANDOM_S_LENS))

.PHONY: test clean random_input

all : nonlazy_protein_benchmark.csv lazy_protein_benchmark.csv

lazy_protein_benchmark.csv :  script_benchmark.py $(FAST_MS)
	python $< --repeat 5 --lazy_wl --output $@ $(BDIR) $(PROTEIN_PREFIXES)

nonlazy_protein_benchmark.csv : script_benchmark.py $(FAST_MS)
	python $< --repeat 5 --output $@ $(BDIR) $(PROTEIN_PREFIXES)

test: script_tests.py $(SLOW_MS) $(FAST_MS)  | $(RANDOM_BDIR)
	python $< --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) $| $(RANDOM_PREFIXES)
	python $< --lazy_wl --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) $| $(RANDOM_PREFIXES)

random_input : generate_input.py | $(RANDOM_BDIR)
	for sl in $(RANDOM_S_LENS); do python $< random $(RANDOM_ALPHABET) --len_t $(RANDOM_T_LEN)  --len_s $$sl --base_dir $|; done

clean  :
	rm -f lazy_protein_benchmark.csv nonlazy_protein_benchmark.csv


