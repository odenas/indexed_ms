
FAST_MS   = fast_ms.py
BDIR      = ./input_data
PREFIXES = $(notdir $(subst .s,,$(wildcard $(BDIR)/*.s)))


stats.csv : $(addsuffix .csv,$(PREFIXES))
	cat $^ | awk -F "," 'NR==1 || $$1!="len_s"' >$@

%.csv : $(FAST_MS) $(BDIR)/%.s $(BDIR)/%.t
	python $^ --load_cst >$@

clean :
	rm -f $(addsuffix .csv,$(PREFIXES))
