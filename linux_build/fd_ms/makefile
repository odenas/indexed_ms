
SDSL_BASE_DIR = /Users/denas/arch/Darwin_x86_64

CXX_WARNINGS = -Wundefined-var-template
INCLUDES = $(SDSL_BASE_DIR)/include
CXXFLAGS = -I$(INCLUDES) -std=c++11 -stdlib=libc++ $(CXX_WARNIGNS) -DNDEBUG -O3

LDFLAGS  = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort64 -ldivsufsort -stdlib=libc++

HEADERS  = cmd_utils.hpp maxrep_vector.hpp stree_sct3.hpp \
		   runs_and_ms_algorithms.hpp utils.hpp input_spec.hpp opt_parser.hpp stree_sct3.hpp\
		   slices.hpp runs_ms.hpp

.PHONY : clean

fd_ms : fd_ms.o
	$(CXX) $(LDFLAGS) $< -o $@


fd_ms.o : fd_ms.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $< -c -o $@


clean : 
	rm -f fd_ms.o fd_ms
