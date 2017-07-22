
SCRIPT = ../../xcode_bin/sandbox_test_maxrep
INPUT  = rep_10Ms_1Mt_abcd


rank_timing.csv : $(addsuffix .csv,$(INPUT))
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.csv : input_data/%.s $(SCRIPT) 
	$(SCRIPT) -load_cst 1 -load_maxrep 1 -s_path $< >$@

clean : 
	rm -f  rank_timing.csv $(addsuffix .csv,$(INPUT))
