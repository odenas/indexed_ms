

FAST_MS          = fast_ms.py
SLOW_MS          = slow_ms.py

ALPHABET         = abcde
RANDOM_BDIR      = datasets/testing
RANDOM_T_LEN     = 200
RANDOM_S_LENS    = 32 64 128 256 512 1024 2048 4096 8192 16384 32768
RANDOM_PREFIXES  = $(addprefix $(ALPHABET)$(RANDOM_T_LEN)_,$(RANDOM_S_LENS))

.PHONY: test_a test_b test_c test_sada clean random_input

all : random_input test_a test_b test_c test_sada

test_a: script_tests.py $(SLOW_MS) $(FAST_MS) | $(RANDOM_BDIR)
	python $< --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) $| $(RANDOM_PREFIXES)

test_b: script_tests.py $(SLOW_MS) $(FAST_MS) | $(RANDOM_BDIR)
	python $< --lazy_wl --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) $| $(RANDOM_PREFIXES)

test_c: script_tests.py $(SLOW_MS) $(FAST_MS) | $(RANDOM_BDIR)
	python $< --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) --nthreads 2 $| $(RANDOM_PREFIXES)

test_sada : script_tests.py $(SLOW_MS) $(FAST_MS) | $(RANDOM_BDIR)
	python $< --sada --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) $| $(RANDOM_PREFIXES)

random_input : generate_input.py $(addsuffix t.txt,$(addprefix $(RANDOM_BDIR)/$(ALPHABET)$(RANDOM_T_LEN)_,$(RANDOM_S_LENS)))

$(RANDOM_BDIR)/$(ALPHABET)$(RANDOM_T_LEN)_%t.txt : generate_input.py | test_input_data
	./generate_input.py random $(ALPHABET) $(subst t.txt,s.txt,$@) $(subst t.txt,,$(subst $(ALPHABET)$(RANDOM_T_LEN)_,,$(notdir $@))) $@  $(RANDOM_T_LEN)


