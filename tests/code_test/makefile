

SLOW_MS     = ./slow_ms.py
BDIR        = ../datasets/testing
RDIR        = ./ms_results
SCRIPT_TEST = ../../py_interface/script_tests.py
PREFIXES    = $(subst .s,,$(notdir $(wildcard $(BDIR)/*.s)))
MS_RESULTS  = $(addsuffix .mstat,$(addprefix $(RDIR)/,$(PREFIXES)))
COMMON_ARGS = --slow_prg $(SLOW_MS) --slow_load_dir $(RDIR) $(BDIR) $(PREFIXES)
#2>/dev/null

.PHONY: test_a test_b test_c test_d test_e

all : test_a test_b test_c1 test_c test_d test_e test_f test_mrep

# lazy=0,fail=0
test_a: $(SCRIPT_TEST) $(SLOW_MS) $(MS_RESULTS)
	python $< $(COMMON_ARGS)

# lazy=1,fail=0
test_b: $(SCRIPT_TEST) $(SLOW_MS)
	python $< --lazy_wl $(COMMON_ARGS) 

# lazy=0,fail=1, maxrep=1
test_c1: $(SCRIPT_TEST) $(SLOW_MS)
	python $< --rank_fail --use_maxrep $(COMMON_ARGS)

# lazy=0,fail=1
test_c: $(SCRIPT_TEST) $(SLOW_MS)
	python $< --rank_fail $(COMMON_ARGS)

# lazy=1,fail=1
test_d: $(SCRIPT_TEST) $(SLOW_MS)
	python $< --lazy_wl --rank_fail $(COMMON_ARGS)

# parallel
test_e: $(SCRIPT_TEST) $(SLOW_MS)
	python $< --nthreads 5 $(COMMON_ARGS)

# lca_parents
test_f : $(SCRIPT_TEST) $(SLOW_MS)
	python $< --lca_parents $(COMMON_ARGS)

# maxrep construction
test_mrep : ../../py_interface/check_maxrep.py 
	python $<  $(BDIR) $(addprefix rnd_200_,32 64 128)

test_speed: $(BDIR)/rnd_100000000s_dis_500000t_abcd.s
	@echo "*** bvectors should take < 150ms ***"
	../fast_ms.py $<  $(subst .s,.t,$<) --load_cst --time_usage \
		| grep bvector

$(MS_RESULTS) : $(RDIR)/%.mstat : $(BDIR)/%.s | $(RDIR)
	python $(SLOW_MS) $< $(subst .s,.t,$<) >$@
