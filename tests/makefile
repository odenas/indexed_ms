

SLOW_MS   = slow_ms.py
BDIR      = datasets/testing
PREFIXES  = $(subst .s,,$(notdir $(wildcard $(BDIR)/*.s)))

.PHONY: test_a test_b test_c

all : test_a test_b test_c

# lazy=0,fail=0
test_a: script_tests.py $(SLOW_MS)
	python $< --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)

# lazy=1,fail=0
test_b: script_tests.py $(SLOW_MS)
	python $< --lazy_wl --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)

# lazy=0,fail=1
test_c: script_tests.py $(SLOW_MS)
	python $< --v --vv --rank_fail --use_maxrep --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)
#python $< --rank_fail --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)

# lazy=1,fail=1
test_d: script_tests.py $(SLOW_MS)
	python $< --lazy_wl --rank_fail --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)

test_e: script_tests.py $(SLOW_MS)
	python $< --slow_prg $(SLOW_MS) --nthreads 5 $(BDIR) $(PREFIXES)


