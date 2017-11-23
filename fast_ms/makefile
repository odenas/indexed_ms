SDSL_BASE_DIR = /home/brt/Documents/projects/indexed_ms/sdsl-lite/build
INCLUDES      = $(SDSL_BASE_DIR)/include	

HEADERS_DIR   = fd_ms
BIN_DIR       = bin

CXXFLAGS      = -I$(INCLUDES) -std=c++11 -ggdb
LDFLAGS       = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lpthread


TARGETS       = $(addprefix $(BIN_DIR)/,$(addsuffix .x,dump_cst dump_maxrep dump_edge_lists dump_nwd_lists sandbox_maxrep sandbox_selectatdist sandbox_parentcalls matching_stats))
OBJECTS       = $(subst .x,.o,$(TARGETS))

all : $(TARGETS)

$(BIN_DIR)/matching_stats.o : matching_stats.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp parent_depth_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_parentcalls.o : sandbox_parentcalls.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp parent_depth_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_maxrep.o : sandbox_maxrep.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp edge_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_selectatdist.o : sandbox_selectatdist.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_nwd_lists.o : dump_nwd_lists.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp parent_depth_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_edge_lists.o : dump_edge_lists.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp maxrep_vector.hpp edge_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_maxrep.o : dump_maxrep.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp maxrep_vector.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_cst.o : dump_cst.cpp  $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@


$(TARGETS) : %.x : %.o
	$(CXX) -o $@ $< $(LDFLAGS)

## clean
clean:
	rm -f $(OBJECTS) $(TARGETS)
