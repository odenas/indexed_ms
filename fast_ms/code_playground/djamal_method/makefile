
SDSL_BASE_DIR = ../../../sdsl-lite/build
INCLUDES      = -I$(SDSL_BASE_DIR)/include	-I../../src

CXXFLAGS      = $(INCLUDES) -std=c++11 -ggdb
LDFLAGS       = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -ldl

TARGETS       = $(addsuffix .x,range_queries)
OBJECTS       = $(subst .x,.o,$(TARGETS))
RLCSA_BIT_O   = $(wildcard ../../src/rlcsa/bits/*.o)
SMSB_BIT_O   = $(addsuffix .o,$(addprefix ../../src/fd_ms/smsb/,range_ms_sum mt19937ar naive_ms_range_sum range_ms_tables))

all : $(TARGETS)

$(OBJECTS) : %.o : %.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

range_queries.x : range_queries.o $(RLCSA_BIT_O) $(SMSB_BIT_O)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean :
	rm -f $(OBJECTS) $(TARGETS) aa.bin
