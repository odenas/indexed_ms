SDSL_BASE_DIR = /home/brt/code/matching_statistics/indexed_ms/sdsl-lite/build

MS_EXEC = ../../../../fast_ms/bin/matching_stats.x
TT_EXEC = ../../../input_stats_data/exploration/index_based_query/node_iterator.x
DUMP_CST= ../../../../fast_ms/bin/dump_ms_data_structures.x

SCRIPT  = mstat_fast_ms.py
OPTS    = --repeat 1 --time_usage --load_cst
ALP     = aaaaaaaaaabbbbbbbbbddddddddddeeeeeeeeeeffffffffffc
#INPUT   = $(foreach m,100 200 400 800 1600 3200,rep_$(m))
INPUT   = $(foreach m,1 500 1000 2000,rep_$(m))
#INPUT   = $(foreach m,2,rep_$(m))

TIME_TARGETS = $(addsuffix .pseq.csv,$(INPUT)) $(addsuffix .lca.csv,$(INPUT))
INPUT_TARGETS = $(addsuffix .s,$(INPUT)) $(addsuffix .s.fwd.stree,$(INPUT)) $(addsuffix .t,$(INPUT))

all : $(TIME_TARGETS)
#all : $(INPUT_TARGETS)

%.lca.csv : %.s.fwd.stree %.t $(MS_EXEC)
	$(SCRIPT) --label lca  --lca_parents $(OPTS) --output $@ $(subst .s.fwd.stree,.s,$<) $(subst .s.fwd.stree,.t,$<)

%.pseq.csv : %.s.fwd.stree %.t $(MS_EXEC)
	$(SCRIPT) --label pseq               $(OPTS) --output $@ $(subst .s.fwd.stree,.s,$<) $(subst .s.fwd.stree,.t,$<)

$(addsuffix .s,$(INPUT)) : %.s :
	mstat_index_input.py --num_blocks 5000 --num_mut $(subst .s,,$(subst rep_,,$@)) rep 10000000 $(ALP) $@

$(addsuffix .t,$(INPUT)) : %.t : %.s $(TT_EXEC)
	./$(TT_EXEC) -load_cst 1 -min_node_depth 10 -max_str_depth 170 -dump_size 100000 -s_path $< >$@

%.s.fwd.stree : %.s $(DUMP_CST)
	$(DUMP_CST) -s_path $< 

clean :
	rm -f stats.csv $(TIME_TARGETS)

clean_data:
	rm -f input.csv parents.csv
	rm -f $(addsuffix .s,$(INPUT)) $(addsuffix .t,$(INPUT))
	rm -f $(addsuffix .s.fwd.stree,$(INPUT)) $(addsuffix .s.rev.stree,$(INPUT)) 
	rm -f  $(addsuffix .s.fwd.maxrep,$(INPUT)) $(addsuffix .s.rev.maxrep,$(INPUT)) 
