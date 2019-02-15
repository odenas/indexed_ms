SDSL_BASE_DIR = /home/odenas/code/matching_statistics/indexed_ms/sdsl-lite/build
INCLUDES      = $(SDSL_BASE_DIR)/include	
INCLUDES_PRIME = /home/cunial/arch/Darwin_x86_64/include

SRC_DIR       = src
HEADERS_DIR   = $(SRC_DIR)/fd_ms
BIN_DIR       = bin
MALLOC_CNT    = malloc_count/malloc_count.o

#CXXFLAGS      = -I$(INCLUDES) -std=c++11 -ggdb
CXXFLAGS      = -I$(INCLUDES)  -I$(INCLUDES_PRIME)    -std=c++11 -O3  -DNDEBUG
LDFLAGS       = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -ldl


TARGET_NAMES  = dump_cst dump_maxrep dump_edge_lists dump_nwd_lists \
				sandbox_maxrep sandbox_selectatdist sandbox_parentcalls \
				input_stats matching_stats \
				matching_stats_parallel freq_profile \
				range_queries range_queries_profile dump_range_index
TARGETS       = $(addprefix $(BIN_DIR)/,$(addsuffix .x,$(TARGET_NAMES)))
OBJECTS       = $(subst .x,.o,$(TARGETS))

all : $(TARGETS)

$(BIN_DIR)/dump_range_index.o : $(SRC_DIR)/dump_range_index.cpp $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter opt_parser input_spec stree_sct3 partial_sums_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/range_queries_profile.o : $(SRC_DIR)/range_queries_profile.cpp $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter opt_parser input_spec stree_sct3 partial_sums_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/range_queries.o : $(SRC_DIR)/range_queries.cpp $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter opt_parser input_spec stree_sct3 partial_sums_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/freq_profile.o : $(SRC_DIR)/freq_profile.cpp $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter opt_parser query))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/matching_stats_parallel.o : $(SRC_DIR)/matching_stats_parallel.cpp $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,slices counter input_spec opt_parser stree_sct3 maxrep_vector query p_runs_vector p_ms_vector ms_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/matching_stats.o : $(SRC_DIR)/matching_stats.cpp $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter input_spec opt_parser stree_sct3 maxrep_vector runs_ms query ms_vector runs_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_parentcalls.o : $(SRC_DIR)/sandbox_parentcalls.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp parent_depth_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_maxrep.o : $(SRC_DIR)/sandbox_maxrep.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp edge_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_selectatdist.o : $(SRC_DIR)/sandbox_selectatdist.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_nwd_lists.o : $(SRC_DIR)/dump_nwd_lists.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp parent_depth_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/input_stats.o : $(SRC_DIR)/input_stats.cpp $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,input_spec opt_parser stree_sct3 stats runs_vector ms_vector maxrep_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_maxrep.o : $(SRC_DIR)/dump_maxrep.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp maxrep_vector.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_edge_lists.o : $(SRC_DIR)/dump_edge_lists.cpp $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp maxrep_vector.hpp edge_list.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_cst.o : $(SRC_DIR)/dump_cst.cpp  $(addprefix $(HEADERS_DIR)/,input_spec.hpp opt_parser.hpp stree_sct3.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@


$(TARGETS) : %.x : %.o $(MALLOC_CNT)
	$(CXX) -o $@ $^ $(LDFLAGS)

## clean
clean:
	rm -f $(OBJECTS) $(TARGETS)
