

BDIR        = ../datasets/testing
RDIR        = ./ms_results
SCRIPT_TEST = ./mstat_script_tests.py
PREFIXES    = $(subst .s,,$(notdir $(wildcard $(BDIR)/rep*.s) $(wildcard $(BDIR)/mut*.s) $(wildcard $(BDIR)/rnd_2*.s)))
MS_RESULTS  = $(addsuffix .mstat,$(addprefix $(RDIR)/,$(PREFIXES)))
COMMON_ARGS = --slow_load_dir $(RDIR) $(BDIR) $(PREFIXES)


.PHONY: test_a test_b test_c test_d test_e test_f

all : test_a test_single_rank test_maxrep test_double_rank

test_a: $(MS_RESULTS)
	$(SCRIPT_TEST) $(COMMON_ARGS)

test_single_rank :
	$(SCRIPT_TEST) $(COMMON_ARGS)
	$(SCRIPT_TEST) $(COMMON_ARGS)               --lazy_wl
	$(SCRIPT_TEST) $(COMMON_ARGS) --lca_parents
	$(SCRIPT_TEST) $(COMMON_ARGS) --lca_parents --lazy_wl

test_maxrep :
	$(SCRIPT_TEST) $(COMMON_ARGS) --lca_parents --rank_fail --use_maxrep_rc      --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS)               --rank_fail --use_maxrep_rc      --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --lca_parents --rank_fail --use_maxrep_vanilla --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS)               --rank_fail --use_maxrep_vanilla --double_rank

test_double_rank :
	$(SCRIPT_TEST) $(COMMON_ARGS)                                     --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS)                         --rank_fail --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS)               --lazy_wl             --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS)               --lazy_wl --rank_fail --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --lca_parents                       --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --lca_parents           --rank_fail --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --lca_parents --lazy_wl             --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --lca_parents --lazy_wl --rank_fail --double_rank


# parallel
test_e:
	@echo "** not working **"
	@echo $(SCRIPT_TEST) --nthreads 5 $(COMMON_ARGS)

test_speed: $(BDIR)/rnd_100000000s_dis_100000t_abcd.s
	@echo "*** bvectors should take < 150ms ***"
	mstat_fast_ms.py $<  $(subst .s,.t,$<) --load_cst --time_usage | grep bvector

$(MS_RESULTS) : $(RDIR)/%.mstat : $(BDIR)/%.s | $(RDIR)
	mstat_slow_ms.py $< $(subst .s,.t,$<) >$@
