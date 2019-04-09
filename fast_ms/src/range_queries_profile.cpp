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
#include "fd_ms/p_ms_vector.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/partial_sums_vector.hpp"
#include "fd_ms/help.hpp"

#include "rlcsa/bits/bitvector.h"
#include "rlcsa/bits/rlevector.h"
#include "rlcsa/bits/deltavector.h"
#include "rlcsa/bits/succinctvector.h"
#include "rlcsa/bits/nibblevector.h"


using namespace std;
using namespace fdms;

typedef StreeOhleb<> cst_t;
typedef typename cst_t::size_type size_type;
typedef typename ms_compression::compression_types Compression;
typedef Counter<size_type> counter_t;

counter_t time_usage{};

class InputFlags {
public:
    size_type block_size;
    size_type range_size, from_idx_max, nqueries;
    bool header;
    Compression compression;

    InputFlags() { }

    InputFlags(const InputFlags& f) :
    block_size{f.block_size},
    range_size{f.range_size}, from_idx_max{f.from_idx_max}, nqueries{f.nqueries},
    header{f.header},
    compression{f.compression} { }

    InputFlags(size_type block_size, size_type range_size, size_type from_idx_max, size_type nqueries, bool header, const Compression compression) :
    block_size{block_size}, range_size{range_size},
    from_idx_max{from_idx_max}, nqueries{nqueries},
    header{header},
    compression{compression} { }

    InputFlags(OptParser input) :
    range_size{static_cast<size_type> (std::stoi(input.getCmdOption("-range_size")))},
    from_idx_max{static_cast<size_type> (std::stoi(input.getCmdOption("-from_max_idx")))},
    nqueries{static_cast<size_type> (std::stoi(input.getCmdOption("-niter")))},
    header{input.getCmdOption("-header") == "1"},
    block_size(static_cast<size_type> (std::stoi(input.getCmdOption("-block_size")))) {
        compression = ms_compression::parse_compression(
            input.getCmdOption("-compression")
        );
    }
};

size_type random_index(size_type max_idx) {
    return static_cast<size_t> (max_idx * static_cast<unsigned long> (std::rand()) / (RAND_MAX + 1UL));
}


template<typename vec_type, typename it_type>
void trivial_comp1(const string ms_path, const InputFlags& flags){
    auto comp_start = timer::now();
    std::ifstream in{ms_path, std::ios::binary};
    vec_type ms(in);
    it_type* it = new it_type(ms);
    time_usage.register_now("load_ms", comp_start);

    time_usage.register_now("load_partial_sums", timer::now());
    comp_start = timer::now();
    for (int k = 0; k < flags.nqueries; k++) {
        size_type start_idx = random_index(flags.from_idx_max);
        size_type end_idx = start_idx + flags.range_size;
        partial_sums_vector1<vec_type, it_type, size_type>::trivial_range_sum(ms, it, start_idx, end_idx);
    }
    time_usage.register_now("algorithm", comp_start);

}

template<typename ms_type, typename ms_sel_1_type>
void trivial_comp(const string ms_path, const InputFlags& flags){
    auto comp_start = timer::now();
    ms_type ms;
    sdsl::load_from_file(ms, ms_path);
    time_usage.register_now("load_ms", comp_start);

    auto ds_start = timer::now();
    ms_sel_1_type ms_sel(&ms);
    time_usage.register_now("select_init", ds_start);


    time_usage.register_now("load_partial_sums", timer::now());

    comp_start = timer::now();
    for (int k = 0; k < flags.nqueries; k++) {
        size_type start_idx = random_index(flags.from_idx_max);
        size_type end_idx = start_idx + flags.range_size;
        partial_sums_vector<size_type, ms_type, ms_sel_1_type>::trivial_range_sum(ms, ms_sel, start_idx, end_idx);
    }
    time_usage.register_now("algorithm", comp_start);
}

template<typename vec_type, typename it_type>
void ridx_comp1(const string ms_path, const string ridx_path, const InputFlags& flags){
    auto comp_start = timer::now();
    std::ifstream in{ms_path, std::ios::binary};
    vec_type ms(in);
    time_usage.register_now("load_ms", comp_start);

    auto ds_start = timer::now();
    it_type* it = new it_type(ms);
    time_usage.register_now("select_init", ds_start);

    ds_start = timer::now();
    sdsl::int_vector<64> ridx;
    sdsl::load_from_file(ridx, ridx_path);
    time_usage.register_now("load_partial_sums", ds_start);

    partial_sums_vector1<vec_type, it_type, size_type> psum(ridx_path, flags.block_size);
    comp_start = timer::now();
    for (int k = 0; k < flags.nqueries; k++) {
        size_type start_idx = random_index(flags.from_idx_max);
        size_type end_idx = start_idx + flags.range_size;
        psum.range_sum(ms, ridx, it, start_idx, end_idx);
    }
    time_usage.register_now("algorithm", comp_start);
}

template<typename ms_type, typename ms_sel_1_type>
void ridx_comp(const string ms_path, const string ridx_path, const InputFlags& flags){
    auto comp_start = timer::now();
    ms_type ms;
    sdsl::load_from_file(ms, ms_path);
    time_usage.register_now("load_ms", comp_start);

    auto ds_start = timer::now();
    ms_sel_1_type ms_sel(&ms);
    time_usage.register_now("select_init", ds_start);

    ds_start = timer::now();
    sdsl::int_vector<64> ridx;
    sdsl::load_from_file(ridx, ridx_path);
    time_usage.register_now("load_partial_sums", ds_start);

    partial_sums_vector<size_type, ms_type, ms_sel_1_type> psum(ridx_path, flags.block_size);
    comp_start = timer::now();
    for (int k = 0; k < flags.nqueries; k++) {
        size_type start_idx = random_index(flags.from_idx_max);
        size_type end_idx = start_idx + flags.range_size;
        psum.range_sum(ms, ridx, ms_sel, start_idx, end_idx);
    }
    time_usage.register_now("algorithm", comp_start);
}

template<typename ms_type, typename ms_sel_1_type>
void comp(const string ms_path, const string ridx_path, const InputFlags& flags) {
    if (flags.block_size > 0) {
        ridx_comp<ms_type, ms_sel_1_type>(ms_path, ridx_path, flags);
    } else {
        trivial_comp<ms_type, ms_sel_1_type>(ms_path, flags);
    }
}

template<typename vec_type, typename enc_type>
void comp1(const string ms_path, const string ridx_path, const InputFlags& flags){
    if (flags.block_size > 0) {
        ridx_comp1<vec_type, enc_type>(ms_path, ridx_path, flags);
    } else {
        trivial_comp1<vec_type, enc_type>(ms_path, flags);
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
    InputFlags flags;
    try{
        flags = InputFlags(input);
    } catch (string s) {
        cerr << s << endl;
        return 1;
    }
    std::srand(123);
    switch(flags.compression)
    {
        case Compression::none:
            comp<sdsl::bit_vector, sdsl::bit_vector::select_1_type>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
            break;
        case Compression::rrr:
            comp<sdsl::rrr_vector<>, sdsl::rrr_vector<>::select_1_type>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
            break;
        case Compression::rle:
            comp1<CSA::RLEVector, CSA::RLEVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
            break;
        case Compression::delta:
            comp1<CSA::DeltaVector, CSA::DeltaVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
            break;
        case Compression::nibble:
            comp1<CSA::NibbleVector, CSA::NibbleVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
            break;
        case Compression::succint:
            comp1<CSA::SuccinctVector, CSA::SuccinctVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
            break;
        default:
            cerr << "Error." << endl;
            break;
    }

    if(flags.header)
        cout << "compression,block_size,range_size,nqueries,method,time_ms" << endl;
    for (auto item : time_usage.reg)
        (cout
            << input.getCmdOption("-compression") << ","
            << flags.block_size << ","
            << flags.range_size << ","
            << flags.nqueries << ","
            << item.first << "," << item.second << endl);
    return 0;
}
