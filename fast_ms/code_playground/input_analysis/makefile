#(.env) brt@brt-ThinkPad-X1-Carbon-4th:buffer_vectors$ g++ -std=c++11 -I ../../sdsl-lite/build/include a.cpp -lsdsl -ldivsufsort  -ldivsufsort64  -L ../../sdsl-lite/build/lib/ -o a.x

SDSL_BASE_DIR = /home/brt/code/matching_statistics/indexed_ms/sdsl-lite/build
INCLUDES      = $(SDSL_BASE_DIR)/include	

DV_SIZE       = 0
CXXFLAGS      = -I$(INCLUDES) -std=c++11 -ggdb
#CXXFLAGS      = -I$(INCLUDES)  -std=c++11 -O3  -DV_SIZE=$(V_SIZE)
LDFLAGS       = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -ldl

TARGETS       = $(addsuffix .x,char_freq)
OBJECTS       = $(subst .x,.o,$(TARGETS))

all : $(TARGETS)

$(OBJECTS) : %.o : %.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(TARGETS) : %.x : %.o 
	$(CXX) -o $@ $^ $(LDFLAGS)

clean :
	rm -f $(OBJECTS) $(TARGETS) aa.bin
