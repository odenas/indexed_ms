

MS_DIR = /home/brt/code/matching_statistics/indexed_ms
IDIR   = ./data
SCRIPT = $(MS_DIR)/fast_ms/bin/matching_stats_parallel.x

INPUT  = $(notdir $(subst .s,,$(wildcard $(IDIR)/protein.s)))
TARGETS= $(foreach nslices,10,$(foreach nthr,2 4 6 8,$(addsuffix .$(nslices).$(nthr).csv,$(INPUT))))
OPTS   = -load_cst 1 -time_usage 1 -nslices 1000


.PHONY: clean

stats.csv: $(TARGETS)
	@echo "plot as a function of slices !!"
	awk -F "," 'BEGIN {print "fname,len_s,len_t,item,value"} {if(FNR > 1) print FILENAME","$$0}' $^ >$@

protein.10.%.csv : $(IDIR)/protein.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads $(subst .csv,,$(subst protein.10.,,$@))  -s_path $< -t_path $(subst .s,,$<).t >$@

protein.100.%.csv : $(IDIR)/protein.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads $(subst .csv,,$(subst protein.100.,,$@))  -s_path $< -t_path $(subst .s,,$<).t >$@

protein.1000.%.csv : $(IDIR)/protein.s $(SCRIPT)
	$(SCRIPT) $(OPTS) -nthreads $(subst .csv,,$(subst protein.1000.,,$@))  -s_path $< -t_path $(subst .s,,$<).t >$@

clean  :
	rm -f stats.csv $(TARGETS)

