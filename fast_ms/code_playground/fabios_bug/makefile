
SDSL_BASE_DIR = ../../../sdsl-lite/build
INCLUDES      = -I$(SDSL_BASE_DIR)/include	-I../../src

CXXFLAGS      = $(INCLUDES) -std=c++11 -ggdb
LDFLAGS       = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -ldl

TARGETS       = $(addsuffix .x,fabio_compress_ms compress_ms)
OBJECTS       = $(subst .x,.o,$(TARGETS))
RLCSA_BIT_O   = $(wildcard ../../src/rlcsa/bits/*.o)
SMSB_BIT_O    = $(addsuffix .o,$(addprefix ../../src/fd_ms/smsb/,range_ms_sum mt19937ar naive_ms_range_sum range_ms_tables))
MALLOC_CNT    = ../../malloc_count/malloc_count.o

all : $(TARGETS)

$(TARGETS) : %.x : %.o $(RLCSA_BIT_O) $(SMSB_BIT_O) $(MALLOC_CNT)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(OBJECTS) : %.o : %.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

clean :
	rm -f $(OBJECTS) $(TARGETS) aa.bin
