
SDSL_BASE_DIR = ../../../../../sdsl-lite/build
SRC_DIR       = ../../../../../fast_ms/src
INCLUDES      = -I$(SDSL_BASE_DIR)/include -I$(SRC_DIR)

CXXFLAGS      = $(INCLUDES) -std=c++11 -O3 -DNDEBUG
#CXXFLAGS      = $(INCLUDES) -std=c++11 -ggdb

LDFLAGS       = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -ldl
MALLOC_CNT    = ../../../../../fast_ms/malloc_count/malloc_count.o ../../../../../fast_ms/malloc_count/stack_count.o

OBJECTS       = rmq_size.o
TARGETS       = rmq_size.x


$(TARGETS) : %.x : %.o $(MALLOC_CNT)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(OBJECTS) : %.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean :
	rm -f $(OBJECTS) $(TARGETS)
