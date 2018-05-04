SDSL_BASE_DIR = /home/brt/code/matching_statistics/indexed_ms/sdsl-lite/build
INCLUDES      = $(SDSL_BASE_DIR)/include	
INCLUDES_PRIME = /home/cunial/arch/Darwin_x86_64/include

SRC_DIR       = src
HEADERS_DIR   = $(SRC_DIR)/fd_ms
BIN_DIR       = bin
MALLOC_CNT    = malloc_count/malloc_count.o

CXXFLAGS      = -I$(INCLUDES) -std=c++11 -ggdb
#CXXFLAGS      = -I$(INCLUDES)  -I$(INCLUDES_PRIME)    -std=c++11 -O3  -DNDEBUG
LDFLAGS       = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -ldl


TARGET_NAMES  = dump_cst dump_maxrep dump_edge_lists dump_nwd_lists \
			sandbox_maxrep sandbox_selectatdist sandbox_parentcalls \
			input_stats matching_stats dump_ms_data_structures runs_tester
TARGETS       = $(addprefix $(BIN_DIR)/,$(addsuffix .x,$(TARGET_NAMES)))
OBJECTS       = $(subst .x,.o,$(TARGETS))

all : $(TARGETS)

$(BIN_DIR)/runs_tester.o : $(SRC_DIR)/runs_tester.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/input_stats.o : $(SRC_DIR)/input_stats.cpp $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,input_spec opt_parser stree_sct3 stats runs_ms maxrep_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/matching_stats.o : $(SRC_DIR)/matching_stats.cpp $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,input_spec opt_parser stree_sct3 maxrep_vector runs_ms query))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_parentcalls.o : $(SRC_DIR)/sandbox_parentcalls.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp parent_depth_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_maxrep.o : $(SRC_DIR)/sandbox_maxrep.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp edge_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_selectatdist.o : $(SRC_DIR)/sandbox_selectatdist.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_nwd_lists.o : $(SRC_DIR)/dump_nwd_lists.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp parent_depth_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_edge_lists.o : $(SRC_DIR)/dump_edge_lists.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp maxrep_vector.hpp edge_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_maxrep.o : $(SRC_DIR)/dump_maxrep.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp maxrep_vector.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_cst.o : $(SRC_DIR)/dump_cst.cpp  $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_ms_data_structures.o : $(SRC_DIR)/dump_ms_data_structures.cpp  $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp maxrep_vector.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

bin/matching_stats.opt.x : $(SRC_DIR)/bin/matching_stats.opt.o
	$(CXX) -o $@ $< $(LDFLAGS)

$(TARGETS) : %.x : %.o $(MALLOC_CNT)
	$(CXX) -o $@ $^ $(LDFLAGS)

## clean
clean:
	rm -f $(OBJECTS) $(TARGETS)
