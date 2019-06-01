

MS_DIR = /home/brt/code/matching_statistics/indexed_ms
IDIR   = ./data
SCRIPT = $(MS_DIR)/fast_ms/bin/matching_stats_parallel.x

INPUT  = $(notdir $(subst .s,,$(wildcard $(IDIR)/*.s)))
TARGETS= $(foreach nslices,10 100 1000,$(foreach nthr,1 2 3 4 5 6 7 8 9 10 11 12 13 14 15,$(addsuffix .$(nslices).$(nthr).csv,$(INPUT))))
OPTS   = -load_cst 1 -time_usage 1 -nslices 1000


.PHONY: clean

stats.csv: $(TARGETS)
	awk -F "," 'BEGIN {print "fname,len_s,len_t,item,value"} {if(FNR > 1) print FILENAME","$$0}' $^ >$@

protein.10.%.csv : $(IDIR)/protein.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads $(subst .csv,,$(subst protein.10.,,$@))  -s_path $< -t_path $(subst .s,,$<).t >$@

protein.100.%.csv : $(IDIR)/protein.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads $(subst .csv,,$(subst protein.100.,,$@))  -s_path $< -t_path $(subst .s,,$<).t >$@

protein.1000.%.csv : $(IDIR)/protein.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads $(subst .csv,,$(subst protein.1000.,,$@))  -s_path $< -t_path $(subst .s,,$<).t >$@

genome.10.%.csv : $(IDIR)/genome.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads $(subst .csv,,$(subst genome.10.,,$@))  -s_path $< -t_path $(subst .s,,$<).t >$@

genome.100.%.csv : $(IDIR)/genome.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads $(subst .csv,,$(subst genome.100.,,$@))  -s_path $< -t_path $(subst .s,,$<).t >$@

genome.1000.%.csv : $(IDIR)/genome.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads $(subst .csv,,$(subst genome.1000.,,$@))  -s_path $< -t_path $(subst .s,,$<).t >$@

clean  :
	rm -f stats.csv $(TARGETS)

