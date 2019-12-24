
SDSL_BASE_DIR = ../../../sdsl-lite/build
SRC_DIR       = ../../src
INCLUDES      = -I$(SDSL_BASE_DIR)/include -I$(SRC_DIR)
CXXFLAGS      = $(INCLUDES) -std=c++11 -ggdb
LDFLAGS       = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -ldl
RLCSA_BIT_O   = $(wildcard $(SRC_DIR)/rlcsa/bits/*.o)
SMSB_BIT_O    = $(wildcard $(SRC_DIR)/fd_ms/smsb/*o)


my_vector.x : my_vector.o $(RLCSA_BIT_O) $(SMSB_BIT_O)
	$(CXX) -o $@ $^ $(LDFLAGS)

my_vector.o : my_vector.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean :
	rm -f $(OBJECTS) $(TARGETS) aa.bin
