

BDIR        = ./input
RDIR        = ./answers

SCRIPT_TEST = ./mstat_script_tests.py
PREFIXES    = $(subst .s,,$(notdir $(wildcard $(BDIR)/*.s)))

# answers obtained with the slow algorithm
MS_RESULTS  = $(addsuffix .mstat,$(addprefix $(RDIR)/,$(PREFIXES)))
# Ms data structures to be used with the fast algorithm
DS_RESULTS  = $(addsuffix .fwd.stree,$(wildcard $(BDIR)/*.s))

COMMON_ARGS = --slow_load_dir $(RDIR) $(BDIR) $(PREFIXES)


.PHONY: test_a test_ds test_single_rank test_maxrep test_double_rank test_e test_range input
all :   test_a test_ds test_single_rank test_maxrep test_double_rank test_e test_range


test_a: input
	$(SCRIPT_TEST) $(COMMON_ARGS) --load_cst --load_maxrep --answer
	$(SCRIPT_TEST) $(COMMON_ARGS) --load_cst --load_maxrep --avg
	$(SCRIPT_TEST) $(COMMON_ARGS) --avg

test_ds : input
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer               --rank_fail --use_maxrep_vanilla --double_rank --load_cst --load_maxrep

test_single_rank : input
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer               --lazy_wl
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer --lca_parents
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer --lca_parents --lazy_wl

test_maxrep : input
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer --lca_parents --rank_fail --use_maxrep_rc      --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer               --rank_fail --use_maxrep_rc      --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer --lca_parents --rank_fail --use_maxrep_vanilla --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer               --rank_fail --use_maxrep_vanilla --double_rank

test_double_rank : input
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer                                     --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer                         --rank_fail --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer               --lazy_wl             --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer               --lazy_wl --rank_fail --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer --lca_parents                       --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer --lca_parents           --rank_fail --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer --lca_parents --lazy_wl             --double_rank
	$(SCRIPT_TEST) $(COMMON_ARGS) --answer --lca_parents --lazy_wl --rank_fail --double_rank

test_range : input
	../bin/dump_range_index.x -ms_path $(BDIR)/rnd_200_128.t.ms  -block_size 8
	for i in 6 7 8 9 10 11 12 13 14; do \
		../bin/range_queries.x  -ms_path $(BDIR)/rnd_200_128.t.ms -ridx_path $(BDIR)/rnd_200_128.t.ms.8.ridx -from_idx 0 -to_idx $$i -block_size 8 -check 1 ;\
	done

test_e: input
	$(SCRIPT_TEST) --nthreads 8 --nslices 4 $(COMMON_ARGS) --answer

input : 
	$(MAKE) -C input

