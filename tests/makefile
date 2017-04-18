

FAST_MS   = fast_ms.py
SLOW_MS   = slow_ms.py
BDIR      = datasets/testing
PREFIXES  = $(subst .s,,$(notdir $(wildcard $(BDIR)/*.s)))

.PHONY: test_a test_b test_c

all : test_a test_b test_c

test_a: script_tests.py $(SLOW_MS) $(FAST_MS)
	python $< --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)

test_b: script_tests.py $(SLOW_MS) $(FAST_MS)
	python $< --lazy_wl --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)

test_c: script_tests.py $(SLOW_MS) $(FAST_MS)
	python $< --fast_prg $(FAST_MS) --slow_prg $(SLOW_MS) --nthreads 5 $(BDIR) $(PREFIXES)


