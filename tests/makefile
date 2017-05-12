

SLOW_MS   = slow_ms.py
BDIR      = datasets/testing
PREFIXES  = $(subst .s,,$(notdir $(wildcard $(BDIR)/*.s)))

.PHONY: test_a test_b test_c test_d test_e

all : test_a test_b test_c1 test_c test_d test_e test_mrep

# lazy=0,fail=0
test_a: script_tests.py $(SLOW_MS)
	python $< --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)

# lazy=1,fail=0
test_b: script_tests.py $(SLOW_MS)
	python $< --lazy_wl --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)

# lazy=0,fail=1
test_c1: script_tests.py $(SLOW_MS)
	python $< --rank_fail --use_maxrep --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)

# lazy=0,fail=1
test_c: script_tests.py $(SLOW_MS)
	python $< --rank_fail --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)

# lazy=1,fail=1
test_d: script_tests.py $(SLOW_MS)
	python $< --lazy_wl --rank_fail --slow_prg $(SLOW_MS) $(BDIR) $(PREFIXES)

test_e: script_tests.py $(SLOW_MS)
	python $< --slow_prg $(SLOW_MS) --nthreads 5 $(BDIR) $(PREFIXES)

test_mrep : check_maxrep.py
	python $<  $(BDIR) $(addprefix rnd_200_,32 64 128)
