NAMES     = mt19937ar naive_ms 
SUM_NAMES = ms_range_sum naive_ms_range_sum ms_range_sum_tables
MAX_NAMES = ms_range_max naive_ms_range_max ms_range_max_tables

TARGETS   = $(addsuffix .o,$(NAMES) $(MAX_NAMES) $(SUM_NAMES))

.PHONY : all
virtual_smsb.a : $(TARGETS)
	ar rcs $@ $^


.PHONY: clean
clean  :
	rm -f $(TARGETS)
