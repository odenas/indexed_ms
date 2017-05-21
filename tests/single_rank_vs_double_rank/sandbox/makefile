
SCRIPT = ../../xcode_bin/single_vs_double_rank
INPUT  = rep_100Ms_500Kt


rank_timing.csv : $(addsuffix .csv,$(INPUT))
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.csv : input_data/%.s $(SCRIPT) 
	$(SCRIPT) -load_cst 1 -s_path $< >$@

clean : 
	rm -f  rank_timing.csv $(addsuffix .csv,$(INPUT))
