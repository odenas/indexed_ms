NAMES     = mt19937ar naive_ms 
SUM_NAMES = ms_range_sum naive_ms_range_sum ms_range_sum_tables
MAX_NAMES = ms_range_max naive_ms_range_max ms_range_max_tables

MAX_SRC     = $(addsuffix .c,$(NAMES) $(MAX_NAMES))
SUM_SRC     = $(addsuffix .c,$(NAMES) $(SUM_NAMES))

MAX_TARGETS = test_ms_range_max test_ms_range_max64
SUM_TARGETS = test_ms_range_sum test_ms_range_sum64


.PHONY : all

all : $(MAX_TARGETS) $(SUM_TARGETS)


$(SUM_TARGETS) : % : %.o $(SUM_SRC)
	gcc $^ -o $@

$(MAX_TARGETS) : % : %.o $(MAX_SRC)
	gcc $^ -o $@


.PHONY: clean
clean  :
	rm -f $(MAX_TARGETS) $(SUM_TARGETS)

