SDSL_BASE_DIR = /home/brt/code/matching_statistics/indexed_ms/sdsl-lite/build

MS_EXEC = ../../../../fast_ms/bin/matching_stats.x
IS_EXEC = ../../../../fast_ms/bin/input_stats.x
TT_EXEC = node_iterator.x

SCRIPT  = mstat_fast_ms.py
OPTS    = --repeat 6 --time_usage --load_cst
ALP     = aaaaaaaaaabbbbbbbbbddddddddddeeeeeeeeeeffffffffffc
#INPUT   = $(foreach m,100 200 400 800 1600 3200,rep_$(m))
INPUT   = $(foreach m,1 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200,rep_$(m))
#INPUT   = $(foreach m,2,rep_$(m))

STAT_TARGETS = $(addsuffix .istat.csv,$(INPUT))
TIME_TARGETS = $(addsuffix .pseq.csv,$(INPUT)) $(addsuffix .lca.csv,$(INPUT))

#all : $(STAT_TARGETS) $(TIME_TARGETS)
#all : $(STAT_TARGETS)
all : $(TIME_TARGETS)

%.istat.csv : %.s %.t
	$(IS_EXEC) -s_path $(subst .istat.csv,.s,$@) -t_path $(subst .istat.csv,.t,$@) -load_cst  1 -load_maxrep 1 >$@

%.lca.csv : %.t $(MS_EXEC)
	$(SCRIPT) --label lca  --lca_parents $(OPTS) --output $@ $(subst .t,.s,$<) $<

%.pseq.csv : %.t $(MS_EXEC)
	$(SCRIPT) --label pseq               $(OPTS) --output $@ $(subst .t,.s,$<) $<

$(addsuffix .s,$(INPUT)) : %.s :
	mstat_index_input.py --num_blocks 5000 --num_mut $(subst .s,,$(subst rep_,,$@)) rep 10000000 $(ALP) $@

$(addsuffix .t,$(INPUT)) : %.t : %.s $(TT_EXEC)
	./$(TT_EXEC) -load_cst 1 -min_node_depth 10 -max_str_depth 170 -dump_size 100000 -s_path $< >$@

$(TT_EXEC) : $(subst .x,.cpp,$(TT_EXEC))
	$(CXX) -I$(SDSL_BASE_DIR)/include -L$(SDSL_BASE_DIR)/lib  -std=c++11 -o $@ $< -lsdsl -ldivsufsort -ldivsufsort64 -ggdb

clean :
	rm -f $(STAT_TARGETS) stats.csv $(TIME_TARGETS)

clean_data:
	rm -f input.csv parents.csv
	rm -f $(addsuffix .s,$(INPUT)) $(addsuffix .t,$(INPUT))
	rm -f $(addsuffix .s.fwd.stree,$(INPUT)) $(addsuffix .s.rev.stree,$(INPUT)) 
	rm -f  $(addsuffix .s.rev.maxrep,$(INPUT)) 
