SDSL_BASE_DIR = ../sdsl-lite/build
INCLUDES      = -I$(SDSL_BASE_DIR)/include

SRC_DIR       = src
HEADERS_DIR   = $(SRC_DIR)/fd_ms
BIN_DIR       = bin
MALLOC_CNT    = malloc_count/malloc_count.o

#CXXFLAGS      = $(INCLUDES) -std=c++11 -O3
CXXFLAGS      = $(INCLUDES)  -std=c++11 -O2  \
				-finline-functions\
				-funswitch-loops\
				-fpredictive-commoning\
				-fgcse-after-reload\
				-ftree-loop-vectorize\
				-ftree-loop-distribute-patterns\
				-fsplit-paths -ftree-slp-vectorize\
				-fvect-cost-model\
				-ftree-partial-pre\
				-fpeel-loops \
				-fipa-cp-clone\
				-DNDEBUG
#-fopt-info-optimized-optall

#CXXFLAGS      = $(INCLUDES) -std=c++11 -ggdb



LDFLAGS       = -L$(SDSL_BASE_DIR)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -ldl
COMMON_DEPS   = $(addprefix $(HEADERS_DIR)/,input_spec.hpp help.hpp opt_parser.hpp)
RLCSA_A       = $(SRC_DIR)/rlcsa/rlcsa.a
DIFF_RLCSA_A  = $(SRC_DIR)/diff/rlcsa/rlcsa.a
LDIFF_RLCSA_A = $(SRC_DIR)/diff_long/rlcsa/rlcsa.a
NDIFF_RLCSA_A = $(SRC_DIR)/diff_none/rlcsa/rlcsa.a
ALL_RLCSA_A   = $(RLCSA_A) $(DIFF_RLCSA_A) $(LDIFF_RLCSA_A) $(NDIFF_RLCSA_A)
SMSB_A        = $(SRC_DIR)/fd_ms/virtual_smsb/virtual_smsb.a

TARGET_NAMES  = dump_cst dump_maxrep dump_edge_lists dump_nwd_lists \
				sandbox_maxrep sandbox_selectatdist sandbox_parentcalls \
				input_stats split_compress_ms \
 				range_queries_profile range_queries \
 				compress_ms diff_compress_ms diff_long_compress_ms diff_none_compress_ms \
 				diff_tables diff_long_tables diff_none_tables \
				matching_stats_parallel matching_stats matching_stats_slow freq_profile \
				dump_range_index print_ms print_int_ms split_ms
TARGETS       = $(addprefix $(BIN_DIR)/,$(addsuffix .x,$(TARGET_NAMES)))
OBJECTS       = $(subst .x,.o,$(TARGETS))

all : $(TARGETS)

# rlcsa and smsb libraries

$(RLCSA_A) : $(wildcard $(dir $(RLCSA_A))/*.cpp) $(wildcard $(dir $(RLCSA_A))/*.h)
	$(MAKE) $(nodir $@) -C $(dir $@)

$(DIFF_RLCSA_A) : $(wildcard $(dir $(DIFF_RLCSA_A))/*.cpp) $(wildcard $(dir $(DIFF_RLCSA_A))/*.h)
	$(MAKE) $(nodir $@) -C $(dir $@)

$(LDIFF_RLCSA_A) : $(wildcard $(dir $(LDIFF_RLCSA_A))/*.cpp) $(wildcard $(dir $(LDIFF_RLCSA_A))/*.h)
	$(MAKE) $(nodir $@) -C $(dir $@)

$(NDIFF_RLCSA_A) : $(wildcard $(dir $(NDIFF_RLCSA_A))/*.cpp) $(wildcard $(dir $(NDIFF_RLCSA_A))/*.h)
	$(MAKE) $(nodir $@) -C $(dir $@)

$(SMSB_A) : $(wildcard $(dir $(SMSB_A))/*.c) $(wildcard $(dir $(SMSB_A))/*.h)
	$(MAKE) $(nodir $@) -C $(dir $@)

# objects

$(BIN_DIR)/matching_stats_slow.o : $(SRC_DIR)/matching_stats_slow.cpp $(COMMON_DEPS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/split_compress_ms.o : $(SRC_DIR)/split_compress_ms.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,slices ms_block_types)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/split_ms.o : $(SRC_DIR)/split_ms.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,slices ms_block_types)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/print_int_ms.o : $(SRC_DIR)/print_int_ms.cpp $(COMMON_DEPS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/print_ms.o : $(SRC_DIR)/print_ms.cpp $(COMMON_DEPS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/compress_ms.o : $(SRC_DIR)/compress_ms.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter p_ms_vector)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/diff_compress_ms.o : $(SRC_DIR)/diff/compress_ms.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter p_ms_vector stree_sct3)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/diff_tables.o : $(SRC_DIR)/diff/createTables.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/diff_long_compress_ms.o : $(SRC_DIR)/diff_long/compress_ms.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter p_ms_vector stree_sct3)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/diff_long_tables.o : $(SRC_DIR)/diff_long/createTables.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/diff_none_compress_ms.o : $(SRC_DIR)/diff_none/compress_ms.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter p_ms_vector stree_sct3)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/diff_none_tables.o : $(SRC_DIR)/diff_none/createTables.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_range_index.o : $(SRC_DIR)/dump_range_index.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 partial_sums_vector partial_max_vector range_query)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/range_queries_profile.o : $(SRC_DIR)/range_queries_profile.cpp $(SRC_DIR)/range_query_commons.hpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter p_ms_vector stree_sct3 partial_sums_vector partial_max_vector range_query)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/range_queries.o : $(SRC_DIR)/range_queries.cpp $(COMMON_DEPS)  $(SRC_DIR)/range_query_commons.hpp $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter p_ms_vector stree_sct3 partial_sums_vector partial_max_vector range_query)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/freq_profile.o : $(SRC_DIR)/freq_profile.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,counter query)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/matching_stats.o : $(SRC_DIR)/matching_stats.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,help slices counter stree_sct3 maxrep_vector query runs_vector ms_vector)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/matching_stats_parallel.o : $(SRC_DIR)/matching_stats_parallel.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,help slices counter stree_sct3 maxrep_vector query p_runs_vector p_ms_vector ms_vector)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_parentcalls.o : $(SRC_DIR)/sandbox_parentcalls.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 parent_depth_list)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_maxrep.o : $(SRC_DIR)/sandbox_maxrep.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 edge_list)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/sandbox_selectatdist.o : $(SRC_DIR)/sandbox_selectatdist.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_nwd_lists.o : $(SRC_DIR)/dump_nwd_lists.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 parent_depth_list)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/input_stats.o : $(SRC_DIR)/input_stats.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 stats runs_vector ms_vector maxrep_vector)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_maxrep.o : $(SRC_DIR)/dump_maxrep.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 maxrep_vector)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_edge_lists.o : $(SRC_DIR)/dump_edge_lists.cpp $(COMMON_DEPS) $(addsuffix .hpp,$(addprefix $(HEADERS_DIR)/,stree_sct3 maxrep_vector edge_list)) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(BIN_DIR)/dump_cst.o : $(SRC_DIR)/dump_cst.cpp $(COMMON_DEPS)  $(addprefix $(HEADERS_DIR)/,stree_sct3.hpp) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -c -o $@

# executables

$(BIN_DIR)/split_compress_ms.x : $(BIN_DIR)/split_compress_ms.o $(RLCSA_A) $(MALLOC_CNT) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)
$(BIN_DIR)/range_queries.x : $(BIN_DIR)/range_queries.o $(RLCSA_A) $(SMSB_A) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)
$(BIN_DIR)/range_queries_profile.x : $(BIN_DIR)/range_queries_profile.o $(RLCSA_A) $(SMSB_A) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/compress_ms.x : $(BIN_DIR)/compress_ms.o $(RLCSA_A) $(MALLOC_CNT) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/diff_compress_ms.x : $(BIN_DIR)/diff_compress_ms.o $(DIFF_RLCSA_A) $(MALLOC_CNT) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/diff_tables.x : $(BIN_DIR)/diff_tables.o | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/diff_long_compress_ms.x : $(BIN_DIR)/diff_long_compress_ms.o $(LDIFF_RLCSA_A) $(MALLOC_CNT) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/diff_long_tables.x : $(BIN_DIR)/diff_long_tables.o | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/diff_none_compress_ms.x : $(BIN_DIR)/diff_none_compress_ms.o $(LDIFF_RLCSA_A) $(MALLOC_CNT) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/diff_none_tables.x : $(BIN_DIR)/diff_none_tables.o | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.x : %.o $(MALLOC_CNT) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR) :
	mkdir $@

## clean
clean:
	rm -f $(OBJECTS) $(TARGETS)

clean-all: clean
	for dr in $(ALL_RLCSA_A) $(SMSB_A); do\
		$(MAKE) clean -C `dirname $$dr`; \
	done

