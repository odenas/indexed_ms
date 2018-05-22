
MS_EXEC = ../../../../fast_ms/bin/matching_stats.x
SCRIPT  = mstat_fast_ms.py
OPTS    = --repeat 6 --time_usage --load_cst
EXEC    = ../../../../fast_ms/bin/input_stats.x
ALP     = aaaaaaaaab \
		  aaaaaaaaabd \
		  aaaaaaaaabcd \
		  aaaaaaaaabcde
SL      = 30000
TL      = 2000
INPUT   = $(foreach a,$(ALP),rep_$(SL)s_dis_$(TL)t_$(a))
TARGETS = $(addsuffix .csv,$(INPUT))
TIME_TARGETS  = $(foreach opts,pseq lca,$(addsuffix .$(opts).csv,$(INPUT)))


#all : echo $(TIME_TARGETS)

rank_timing.csv : $(TIME_TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

stats.csv : $(TARGETS)
	awk -F "," 'NR==1 || $$1!="len_s"{print $$0","FILENAME}' $^ >$@

$(TARGETS) : %.csv : %.s $(EXEC) 
	$(EXEC) -s_path $<  -t_path $(subst .s,.t,$<) -load_maxrep 0 -load_cst 1 >$@

$(addsuffix .s,$(INPUT)) : %.s : generate_input3.py
	python $< rep_dis $(subst .s,,$(subst rep_$(SL)s_dis_$(TL)t_,,$@)) $(SL) $(TL)  --rep_blocks 100 --rep_mut 30

%.lca.csv : %.s $(MS_EXEC)
	$(SCRIPT) --label lca  --lca_parents $(OPTS) --output $@ $< $(subst .s,,$<).t

%.pseq.csv : %.s $(MS_EXEC)
	$(SCRIPT) --label pseq               $(OPTS) --output $@ $< $(subst .s,,$<).t

clean :
	rm -f $(TARGETS) stats.csv $(TIME_TARGETS)

clean_data:
	rm -f $(addsuffix .s,$(INPUT)) $(addsuffix .t,$(INPUT)) \
		$(addsuffix .s.fwd.stree,$(INPUT)) $(addsuffix .s.rev.stree,$(INPUT))
