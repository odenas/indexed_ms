
FAST_MS        = ../fast_ms.py
BDIR           = input_data
SINGLE_TARGETS = $(addsuffix .csv,l0_f0_m0 l0_f0_m1 l0_f1_m0 l0_f1_m1 l1_f0_m0 l1_f0_m1 l1_f1_m0 l1_f1_m1)
INPUT          = $(addprefix $(BDIR)/rnd_1Gs_500Kt.,s t)
COMMON_ARGS     = --space_usage --time_usage --load_cst



res.csv : $(SINGLE_TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

# lazy=0,fail=0,maxrep0
l0_f0_m0.csv: $(FAST_MS)
	python $< $(COMMON_ARGS) $(INPUT)\
		| awk -F "," '{if(NR == 1) {print $$0",label"} else {print $$0",$(subst .csv,,$@)"}}' >$@

l0_f0_m1.csv: $(FAST_MS)
	python $< $(COMMON_ARGS)  --use_maxrep $(INPUT)\
		| awk -F "," '{if(NR == 1) {print $$0",label"} else {print $$0",$(subst .csv,,$@)"}}' >$@

l0_f1_m0.csv: $(FAST_MS)
	python $< $(COMMON_ARGS)  --rank_fail $(INPUT)\
		| awk -F "," '{if(NR == 1) {print $$0",label"} else {print $$0",$(subst .csv,,$@)"}}' >$@

l0_f1_m1.csv: $(FAST_MS)
	python $< $(COMMON_ARGS)  --use_maxrep --rank_fail $(INPUT)\
		| awk -F "," '{if(NR == 1) {print $$0",label"} else {print $$0",$(subst .csv,,$@)"}}' >$@

l1_f0_m0.csv: $(FAST_MS)
	python $< $(COMMON_ARGS)  --lazy_wl $(INPUT)\
		| awk -F "," '{if(NR == 1) {print $$0",label"} else {print $$0",$(subst .csv,,$@)"}}' >$@

l1_f0_m1.csv: $(FAST_MS)
	python $< $(COMMON_ARGS)  --use_maxrep --lazy_wl $(INPUT)\
		| awk -F "," '{if(NR == 1) {print $$0",label"} else {print $$0",$(subst .csv,,$@)"}}' >$@

l1_f1_m0.csv: $(FAST_MS)
	python $< $(COMMON_ARGS)  --lazy_wl --rank_fail $(INPUT)\
		| awk -F "," '{if(NR == 1) {print $$0",label"} else {print $$0",$(subst .csv,,$@)"}}' >$@

l1_f1_m1.csv: $(FAST_MS)
	python $< $(COMMON_ARGS)  --use_maxrep --lazy_wl --rank_fail $(INPUT)\
		| awk -F "," '{if(NR == 1) {print $$0",label"} else {print $$0",$(subst .csv,,$@)"}}' >$@

clean :
	rm -f $(SINGLE_TARGETS) res.csv


