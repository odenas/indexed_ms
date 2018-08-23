

MS_DIR = /home/brt/code/matching_statistics/indexed_ms
IDIR   = $(MS_DIR)/experiments/datasets/big_paper3
SCRIPT = $(MS_DIR)/fast_ms/bin/matching_stats_parallel.x

INPUT  = $(notdir $(subst .s,,$(wildcard $(IDIR)/rnd_100000000s_dis_500000t_abcde*.s) $(wildcard $(IDIR)/rnd_100000000s_sim_500000t_abcde*.s)))
TARGETS= $(foreach opts,1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32,$(addsuffix .$(opts).csv,$(INPUT)))
OPTS   = -load_cst 1 -time_usage 1


.PHONY: clean

stats.csv: $(TARGETS)
	awk -F "," 'BEGIN {print "fname,len_s,len_t,item,value"} {if(FNR > 1) print FILENAME","$$0}' $^ >$@

%.1.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 1  -s_path $< -t_path $(subst .s,,$<).t >$@

%.2.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 2  -s_path $< -t_path $(subst .s,,$<).t >$@

%.3.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 3  -s_path $< -t_path $(subst .s,,$<).t >$@

%.4.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 4  -s_path $< -t_path $(subst .s,,$<).t >$@

%.5.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 5  -s_path $< -t_path $(subst .s,,$<).t >$@

%.6.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 6  -s_path $< -t_path $(subst .s,,$<).t >$@

%.7.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 7  -s_path $< -t_path $(subst .s,,$<).t >$@

%.8.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 8  -s_path $< -t_path $(subst .s,,$<).t >$@

%.9.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 9  -s_path $< -t_path $(subst .s,,$<).t >$@

%.10.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 10  -s_path $< -t_path $(subst .s,,$<).t >$@

%.11.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 11  -s_path $< -t_path $(subst .s,,$<).t >$@

%.12.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 12  -s_path $< -t_path $(subst .s,,$<).t >$@

%.13.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 13  -s_path $< -t_path $(subst .s,,$<).t >$@

%.14.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 14  -s_path $< -t_path $(subst .s,,$<).t >$@

%.15.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 15  -s_path $< -t_path $(subst .s,,$<).t >$@

%.16.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 16 -s_path $< -t_path $(subst .s,,$<).t >$@

%.17.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 17 -s_path $< -t_path $(subst .s,,$<).t >$@

%.18.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 18 -s_path $< -t_path $(subst .s,,$<).t >$@

%.19.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 19 -s_path $< -t_path $(subst .s,,$<).t >$@

%.20.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 20 -s_path $< -t_path $(subst .s,,$<).t >$@

%.21.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 21 -s_path $< -t_path $(subst .s,,$<).t >$@

%.22.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 22 -s_path $< -t_path $(subst .s,,$<).t >$@

%.23.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 23 -s_path $< -t_path $(subst .s,,$<).t >$@

%.24.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 24 -s_path $< -t_path $(subst .s,,$<).t >$@

%.25.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 25 -s_path $< -t_path $(subst .s,,$<).t >$@

%.26.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 26 -s_path $< -t_path $(subst .s,,$<).t >$@

%.27.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 27 -s_path $< -t_path $(subst .s,,$<).t >$@

%.28.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 28 -s_path $< -t_path $(subst .s,,$<).t >$@

%.29.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 29 -s_path $< -t_path $(subst .s,,$<).t >$@

%.30.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 30 -s_path $< -t_path $(subst .s,,$<).t >$@

%.31.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 31 -s_path $< -t_path $(subst .s,,$<).t >$@

%.32.csv : $(IDIR)/%.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads 32 -s_path $< -t_path $(subst .s,,$<).t >$@

clean  :
	rm -f stats.csv $(TARGETS)
