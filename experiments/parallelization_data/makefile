

MS_DIR = /home/brt/code/matching_statistics/indexed_ms
IDIR   = $(MS_DIR)/experiments/datasets/big_paper3
SCRIPT = $(MS_DIR)/fast_ms/bin/matching_stats_parallel.x

INPUT  = $(notdir $(subst .s,,$(wildcard $(IDIR)/rnd_100000000s_dis_500000t_abcde*.s) $(wildcard $(IDIR)/rnd_100000000s_sim_500000t_abcde*.s)))
TARGETS= $(foreach opts,1 2 4 8 16 32,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS   = -load_cst 1 -time_usage 1


.PHONY: clean

stats.csv: $(TARGETS)
	awk -F "," 'BEGIN {print "fname,len_s,len_t,item,value"} {if(FNR > 1) print FILENAME","$$0}' $^ >$@

%.1.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 1  -s_path $< -t_path $(subst .s,,$<).t >$@

%.2.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 2  -s_path $< -t_path $(subst .s,,$<).t >$@

%.4.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 4  -s_path $< -t_path $(subst .s,,$<).t >$@

%.8.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 8  -s_path $< -t_path $(subst .s,,$<).t >$@

%.16.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 16 -s_path $< -t_path $(subst .s,,$<).t >$@

%.32.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 32 -s_path $< -t_path $(subst .s,,$<).t >$@

clean  :
	rm -f stats.csv $(TARGETS)
