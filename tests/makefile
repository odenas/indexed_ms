
include makefile.paths

ALPHABET         = abcde
RANDOM_BDIR      = test_input_data
RANDOM_T_LEN     = 200
RANDOM_S_LENS    = 32 64 128 256 512 1024 2048 4096 8192 16384 32768
RANDOM_PREFIXES  = $(addprefix $(ALPHABET)$(RANDOM_T_LEN)_,$(RANDOM_S_LENS))

.PHONY: test_a test_b clean random_input

all : random_input  test_a test_b

test_a: script_tests.py $(SLOW_MS) $(FAST_MS) | $(RANDOM_BDIR)
	python $< --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) $| $(RANDOM_PREFIXES)

test_b: script_tests.py $(SLOW_MS) $(FAST_MS) | $(RANDOM_BDIR)
	python $< --lazy_wl --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) $| $(RANDOM_PREFIXES)

random_input : generate_input.py $(addsuffix t.txt,$(addprefix $(RANDOM_BDIR)/$(ALPHABET)$(RANDOM_T_LEN)_,$(RANDOM_S_LENS)))

$(RANDOM_BDIR)/$(ALPHABET)$(RANDOM_T_LEN)_%t.txt : generate_input.py | test_input_data
	./generate_input.py random $(ALPHABET) --base_dir $(RANDOM_BDIR) --len_t $(RANDOM_T_LEN) --len_s $(subst t.txt,,$(subst $(ALPHABET)$(RANDOM_T_LEN)_,,$(notdir $@)))


clean  :
	rm -f $(RANDOM_BDIR)/*txt


