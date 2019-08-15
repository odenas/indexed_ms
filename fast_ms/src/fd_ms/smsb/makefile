NAMES   = range_ms_sum\
		  mt19937ar\
		  naive_ms_range_sum\
		  range_ms_tables

SRC     = $(addsuffix .c,$(NAMES))
TARGETS = test_range_ms64 test_range_ms


.PHONY : all

all : $(TARGETS)


$(TARGETS) : % : %.c $(SRC)
	gcc $^ -o $@


.PHONY: clean
clean  :
	rm -f $(TARGETS)
