
FAST_MS   = fast_ms.py
PREFIXES = mut_100Ms_5Mt_10 mut_100Ms_5Mt_100 mut_100Ms_5Mt_1000 mut_100Ms_5Mt_10000 rnd_100Ms_5Mt rep_100Ms_5Mt


stats.csv : $(addsuffix .csv,$(PREFIXES))
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.csv : $(FAST_MS) input_data/%.s input_data/%.t
	python $^ --load_maxrep --load_cst >$@

clean :
	rm -f $(addsuffix .csv,$(PREFIXES))
