
SDSL_BASE_DIR = ../../../sdsl-lite/build
INCLUDES      = $(SDSL_BASE_DIR)/include

SRC_DIR       = ../../../fast_ms/src
HEADERS_DIR   = $(SRC_DIR)/fd_ms
BIN_DIR       = bin
MALLOC_CNT    = $(SRC_DIR)/../malloc_count/malloc_count.o

#CXXFLAGS      = -I$(INCLUDES) -std=c++11 -ggdb
CXXFLAGS      = -I$(INCLUDES)  -I$(SRC_DIR) -I$(SRC_DIR)/.. -std=c++11 -O3  -DNDEBUG -ffast-math -funroll-loops
LDFLAGS       = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -ldl
COMMON_DEPS   = $(addprefix $(HEADERS_DIR)/,input_spec.hpp help.hpp opt_parser.hpp)
RLCSA_BIT_O   = $(wildcard $(SRC_DIR)/rlcsa/bits/*.o)

TARGETS  = memory_usage.x
OBJECTS  = $(subst .x,.o,$(TARGETS))

all : $(TARGETS)

memory_usage.o : memory_usage.cpp $(COMMON_DEPS)
	$(CXX) $(CXXFLAGS) $< -c -o $@

# executables

memory_usage.x : memory_usage.o $(RLCSA_BIT_O) $(MALLOC_CNT)
	$(CXX) -o $@ $^ $(LDFLAGS)

## clean
clean:
	rm -f $(OBJECTS) $(TARGETS)

