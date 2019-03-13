#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <future>
#include <thread>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/select_support.hpp>

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/counter.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/partial_sums_vector.hpp"
#include "fd_ms/help.hpp"

using namespace std;
using namespace fdms;

typedef typename StreeOhleb<>::size_type size_type;
#ifdef COMPRESSED
typedef typename sdsl::rrr_vector<> ms_type;
typedef typename sdsl::rrr_vector<>::select_1_type ms_sel_1_type;
#else
typedef typename sdsl::bit_vector ms_type;
typedef typename sdsl::bit_vector::select_1_type ms_sel_1_type;
#endif


typedef Counter<size_type> counter_t;

counter_t time_usage{};

class InputFlags {
public:
    size_type block_size;
    size_type range_size, from_idx_max, nqueries;
    bool header;

    InputFlags() { }

    InputFlags(const InputFlags& f) :
    block_size{f.block_size},
    range_size{f.range_size}, from_idx_max{f.from_idx_max}, nqueries{f.nqueries},
    header{f.header} { }

    InputFlags(size_type block_size, size_type range_size, size_type from_idx_max, size_type nqueries, bool header) :
    block_size{block_size}, range_size{range_size},
    from_idx_max{from_idx_max}, nqueries{nqueries},
    header{header} { }

    InputFlags(OptParser input) :
    range_size{static_cast<size_type> (std::stoi(input.getCmdOption("-range_size")))},
    from_idx_max{static_cast<size_type> (std::stoi(input.getCmdOption("-from_max_idx")))},
    nqueries{static_cast<size_type> (std::stoi(input.getCmdOption("-niter")))},
    header{input.getCmdOption("-header") == "1"},
    block_size(static_cast<size_type> (std::stoi(input.getCmdOption("-block_size")))) {
    }
};

size_type random_index(size_type max_idx) {
    return static_cast<size_t> (max_idx * static_cast<unsigned long> (std::rand()) / (RAND_MAX + 1UL));
}

void trivial_comp(ms_type& ms, ms_sel_1_type& ms_sel, const InputFlags& flags){

    time_usage.register_now("load_partial_sums", timer::now());

    auto comp_start = timer::now();
    for (int k = 0; k < flags.nqueries; k++) {
        size_type start_idx = random_index(flags.from_idx_max);
        size_type end_idx = start_idx + flags.range_size;
        partial_sums_vector<size_type, ms_type, ms_sel_1_type>::trivial_range_sum(ms, ms_sel, start_idx, end_idx);
    }
    time_usage.register_now("algorithm", comp_start);
}

void ridx_comp(ms_type& ms, ms_sel_1_type& ms_sel, const string ridx_path, const InputFlags& flags){
    auto ds_start = timer::now();
    sdsl::int_vector<64> ridx;
    sdsl::load_from_file(ridx, ridx_path);
    time_usage.register_now("load_partial_sums", ds_start);

    partial_sums_vector<size_type, ms_type, ms_sel_1_type> psum(ridx_path, flags.block_size);
    auto comp_start = timer::now();
    for (int k = 0; k < flags.nqueries; k++) {
        size_type start_idx = random_index(flags.from_idx_max);
        size_type end_idx = start_idx + flags.range_size;
        psum.range_sum(ms, ridx, ms_sel, start_idx, end_idx);
    }
    time_usage.register_now("algorithm", comp_start);
}

void comp(const string ms_path, const string ridx_path, const InputFlags& flags) {
    auto comp_start = timer::now();
    ms_type ms;
    sdsl::load_from_file(ms, ms_path);
    time_usage.register_now("load_ms", comp_start);

    auto ds_start = timer::now();
    ms_sel_1_type ms_sel(&ms);
    time_usage.register_now("select_init", ds_start);

    std::srand(123);

    size_type answer = 0;
    if (flags.block_size > 0) {
        ridx_comp(ms, ms_sel, ridx_path, flags);
    } else {
        trivial_comp(ms, ms_sel, flags);
    }
}

int main(int argc, char **argv) {
    if(argc == 1){
        (cerr << "Run a series of range queries with random end points and report (in csv format) their time in the standard output.\n"
              << "Args:\n"
              << help__ms_path
              << help__ridx_path
              << help_block_size
              << "\t-range_size <positive int>: the size of the range queries\n"
              << "\t-from_max_idx <non-negative int>: range starting point will be greater or equal to this value\n"
              << "\t-niter <positive int>: the number of queries to run\n"
              << "\t-header 1: print a header of the report\n"
              << endl);
        exit(0);
    }
    OptParser input(argc, argv);
    InputFlags flags(input);
    comp(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);

    if(flags.header)
        cout << "block_size,range_size,nqueries,method,time_ms" << endl;
    for (auto item : time_usage.reg)
        (cout << flags.block_size << "," << flags.range_size << "," << flags.nqueries << ","
            << item.first << "," << item.second << endl);
    return 0;
}
