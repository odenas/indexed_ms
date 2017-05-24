
SCRIPT = script_benchmark.py
INPUT  = mut_100Ms_5Mt_10 mut_100Ms_5Mt_100  \
		 mut_100Ms_5Mt_1000 mut_100Ms_5Mt_10000 mut_100Ms_5Mt_100000 \
		 mut_100Ms_5Mt_20 mut_100Ms_5Mt_30 mut_100Ms_5Mt_40 mut_100Ms_5Mt_50 \
		 mut_100Ms_5Mt_500000 mut_100Ms_5Mt_60 mut_100Ms_5Mt_70 \
		 mut_100Ms_5Mt_80 mut_100Ms_5Mt_90 rnd_100Ms_5Mt rep_100Ms_5Mt
REPEAT = 8


res.csv : double_rank.csv double_rank_maxrep.csv
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

double_rank.csv :  ./script_benchmark.py $(addsuffix .s,$(addprefix input_data/,$(INPUT)))
	python $< --load_cst --repeat $(REPEAT) --output $@ --label $(subst .csv,,$@) input_data $(INPUT)

double_rank_maxrep.csv :  ./script_benchmark.py $(addsuffix .s,$(addprefix input_data/,$(INPUT)))
	python $< --load_cst --use_maxrep --load_maxrep --repeat $(REPEAT) --output $@ --label $(subst .csv,,$@) input_data $(INPUT)

clean  : 
	rm -f res.csv double_rank_maxrep.csv double_rank.csv
