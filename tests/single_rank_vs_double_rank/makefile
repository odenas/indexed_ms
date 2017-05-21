
SCRIPT = script_benchmark.py
INPUT  = mut_100Ms_5Mt_10 mut_100Ms_5Mt_100  \
		 mut_100Ms_5Mt_1000 mut_100Ms_5Mt_10000 mut_100Ms_5Mt_100000 \
		 mut_100Ms_5Mt_20 mut_100Ms_5Mt_30 mut_100Ms_5Mt_40 mut_100Ms_5Mt_50 \
		 mut_100Ms_5Mt_500000 mut_100Ms_5Mt_60 mut_100Ms_5Mt_70 \
		 mut_100Ms_5Mt_80 mut_100Ms_5Mt_90 rnd_100Ms_5Mt rep_100Ms_5Mt


rank_timing.csv : $(addsuffix .csv,$(INPUT)) $(addsuffix .fail.csv,$(INPUT))
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.fail.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --label rank_fail --rank_fail --load_cst --repeat 5 input_data $(subst .fail.csv,,$@) >$@

%.csv : input_data/%.s $(SCRIPT) 
	python $(SCRIPT) --load_cst --repeat 5 input_data $(subst .csv,,$@) >$@

clean : 
	rm -f  rank_timing.csv $(addsuffix .csv,$(INPUT)) $(addsuffix .fail.csv,$(INPUT))
