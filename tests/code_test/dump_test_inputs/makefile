
DUMP_CST      = ../../../fast_ms/bin/dump_cst.x
DUMP_MAXREP   = ../../../fast_ms/bin/dump_maxrep.x
M_STAT        =  ../../../fast_ms/bin/matching_stats.x
INPUTS        = $(subst .s,,$(wildcard *.s))

all : $(addsuffix .ms0,$(INPUTS)) $(addsuffix .ms1,$(INPUTS))
	for f in $(INPUTS); \
	do \
		diff -q $$f.ms0 $$f.ms1; \
	done


%.ms1 : %.s %.t %.fwd.stree %.rev.maxrep
	$(M_STAT) \
		-s_path $(subst ms1,s,$@) \
		-t_path $(subst ms1,t,$@) \
		-load_cst 1 -load_maxrep 1 \
		-rank_fail 1 -double_rank 1 -use_maxrep_rc 1 \
		-answer 1 >$@

%.ms0 : %.s %.t
	$(M_STAT) \
		-s_path $(subst ms0,s,$@) \
		-t_path $(subst ms0,t,$@) \
		-rank_fail 1 -double_rank 1 -use_maxrep_rc 1 \
		-answer 1 >$@

%.fwd.stree : %.s $(DUMP_CST)
	$(DUMP_CST) -s_path $< 

%.rev.maxrep : %.s %.fwd.stree $(DUMP_MAXREP)
	$(DUMP_MAXREP) -load_cst 1 -s_path $<

clean :
	rm -f $(addsuffix .ms0,$(INPUTS)) \
		$(addsuffix .ms1,$(INPUTS)) \
		$(addsuffix .rev.maxrep,$(INPUTS)) \
		$(addsuffix .rev.stree,$(INPUTS)) \
		$(addsuffix .fwd.stree,$(INPUTS))
