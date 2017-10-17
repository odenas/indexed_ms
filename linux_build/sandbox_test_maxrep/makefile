SDSL_BASE_DIR = /Users/denas/arch/Darwin_x86_64

CXX_WARNINGS = -Wundefined-var-template
INCLUDES = $(SDSL_BASE_DIR)/include
CXXFLAGS = -I$(INCLUDES) -std=c++11 -stdlib=libc++ $(CXX_WARNIGNS) -DNDEBUG -O3

LDFLAGS  = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort64 -ldivsufsort -stdlib=libc++

HEADERS  = cmd_utils.hpp maxrep_construction.hpp stree_sct3.hpp runs_and_ms_algorithms.hpp utils.hpp cst_iterator.hpp

.PHONY : clean

main : main.o
	$(CXX) $(LDFLAGS) $< -o $@

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $< -c -o $@


clean : 
	rm -f main.o main

