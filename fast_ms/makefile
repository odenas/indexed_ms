SDSL_BASE_DIR = ../sdsl-lite/build
INCLUDES      = $(SDSL_BASE_DIR)/include	

SRC_DIR       = src
HEADERS_DIR   = $(SRC_DIR)/fd_ms
BIN_DIR       = bin
MALLOC_CNT    = malloc_count/malloc_count.o

#CXXFLAGS      = -I$(INCLUDES) -std=c++11 -ggdb
CXXFLAGS      = -I$(INCLUDES)  -std=c++11 -O3  -DNDEBUG
LDFLAGS       = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -ldl
COMMON_DEPS   = $(addprefix $(HEADERS_DIR)/,input_spec.hpp help.hpp opt_parser.hpp)
RLCSA_BIT_O   = $(wildcard $(SRC_DIR)/rlcsa/bits/*.o)


TARGET_NAMES  = dump_cst dump_maxrep dump_edge_lists dump_nwd_lists \
				sandbox_maxrep sandbox_selectatdist sandbox_parentcalls \
				input_stats matching_stats \
				matching_stats_parallel freq_profile \
				dump_range_index
TARGETS       = $(addprefix $(BIN_DIR)/,$(addsuffix .x,$(TARGET_NAMES)))
OBJECTS       = $(subst .x,.o,$(TARGETS))

all : $(TARGETS) $(addprefix $(BIN_DIR)/,compress_ms.x range_queries_profile.x range_queries.x)

$(BIN_DIR)/compress_ms.o : $(SRC_DIR)/compress_ms.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter p_ms_vector partial_sums_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_range_index.o : $(SRC_DIR)/dump_range_index.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter stree_sct3 partial_sums_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/range_queries_profile.o : $(SRC_DIR)/range_queries_profile.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter p_ms_vector stree_sct3 partial_sums_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/range_queries.o : $(SRC_DIR)/range_queries.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter p_ms_vector stree_sct3 partial_sums_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/freq_profile.o : $(SRC_DIR)/freq_profile.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter query))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/matching_stats_parallel.o : $(SRC_DIR)/matching_stats_parallel.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,slices counter stree_sct3 maxrep_vector query p_runs_vector p_ms_vector ms_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/matching_stats.o : $(SRC_DIR)/matching_stats.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter stree_sct3 maxrep_vector runs_ms query ms_vector runs_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_parentcalls.o : $(SRC_DIR)/sandbox_parentcalls.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 parent_depth_list))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_maxrep.o : $(SRC_DIR)/sandbox_maxrep.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 edge_list))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_selectatdist.o : $(SRC_DIR)/sandbox_selectatdist.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_nwd_lists.o : $(SRC_DIR)/dump_nwd_lists.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 parent_depth_list))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/input_stats.o : $(SRC_DIR)/input_stats.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 stats runs_vector ms_vector maxrep_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_maxrep.o : $(SRC_DIR)/dump_maxrep.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 maxrep_vector))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_edge_lists.o : $(SRC_DIR)/dump_edge_lists.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 maxrep_vector edge_list))
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_cst.o : $(SRC_DIR)/dump_cst.cpp $(COMMON_DEPS)  $(addprefix $(HEADERS_DIR)/,stree_sct3.hpp)
	$(CXX) $(CXXFLAGS) $< -c -o $@

# executables

$(BIN_DIR)/range_queries.x : $(BIN_DIR)/range_queries.o $(RLCSA_BIT_O)
	$(CXX) -o $@ $^ $(LDFLAGS)
$(BIN_DIR)/range_queries_profile.x : $(BIN_DIR)/range_queries_profile.o $(RLCSA_BIT_O)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/compress_ms.x : $(BIN_DIR)/compress_ms.o $(RLCSA_BIT_O)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(TARGETS) : %.x : %.o $(MALLOC_CNT)
	$(CXX) -o $@ $^ $(LDFLAGS)

## clean
clean:
	rm -f $(OBJECTS) $(TARGETS)
