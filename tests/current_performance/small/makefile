
FAST_MS     = ../../fast_ms.py
BDIR        = input_data
INPUT       = $(addprefix $(BDIR)/rep_1000000s_sim_200000t_abcdefghijklmnopqrst_sim1000.,s t)
TARGETS     = $(addsuffix .csv,l0_f0_m0 l0_f1_m0 l0_f1_m1 l1_f0_m0 l1_f1_m0)
COMMON_ARGS = --space_usage --time_usage --load_cst


res.csv : $(TARGETS)
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

# lazy=0,fail=0,maxrep0
l0_f0_m0.csv: $(FAST_MS) $(INPUT)
	python $^ $(COMMON_ARGS) --label $(subst .csv,,$@) --output $@

#%.l0_f0_m1.csv : NA

l0_f1_m0.csv: $(FAST_MS) $(INPUT)
	python $^ $(COMMON_ARGS) --label $(subst .csv,,$@) --rank_fail --output $@

l0_f1_m1.csv: $(FAST_MS) $(INPUT)
	python $^ $(COMMON_ARGS) --label $(subst .csv,,$@) --rank_fail --use_maxrep --output $@

l1_f0_m0.csv: $(FAST_MS) $(INPUT)
	python $^ $(COMMON_ARGS) --label $(subst .csv,,$@) --lazy_wl --output $@

# l1_f0_m1.csv: NA 

l1_f1_m0.csv: $(FAST_MS) $(INPUT)
	python $^ $(COMMON_ARGS) --label $(subst .csv,,$@) --rank_fail --lazy_wl --output $@

#l1_f1_m1.csv: NA

clean :
	rm -f $(TARGETS) res.csv
