NAMES = range_ms_sum gen_range_sum_lookup_tables  mt19937ar  naive_ms_range_sum  range_ms_tables test_range_ms64 test_range_ms
OBJ = $(addsuffix .o,$(NAMES))

all: $(OBJ)

$(OBJECTS) : %.o : %.c
	$(CXX) -Wall $< -c -o $@
