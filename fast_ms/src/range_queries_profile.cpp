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
#include "range_query_commons.hpp"

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


class InputFlags {
public:
    int64_t block_size;
    size_type range_size, from_idx_max, nqueries;
    bool header;
    Compression compression;

    InputFlags() { }

    InputFlags(const InputFlags& f) :
    block_size{f.block_size},
    range_size{f.range_size}, from_idx_max{f.from_idx_max}, nqueries{f.nqueries},
    header{f.header},
    compression{f.compression} { }

    InputFlags(OptParser input) :
    range_size{static_cast<size_type> (std::stoi(input.getCmdOption("-range_size")))},
    from_idx_max{static_cast<size_type> (std::stoi(input.getCmdOption("-from_max_idx")))},
    nqueries{static_cast<size_type> (std::stoi(input.getCmdOption("-niter")))},
    header{input.getCmdOption("-header") == "1"},
    block_size(static_cast<int64_t> (std::stoi(input.getCmdOption("-block_size")))) {
        compression = ms_compression::parse_compression(
            input.getCmdOption("-compression")
        );
    }
};

template<typename dispatcher_t>
void rle_comp(const string& ms_path, const string& ridx_path, rq_dispatcher::counter_t& time_usage, const InputFlags& flags){
    if(flags.block_size == 0)
        return dispatcher_t::trivial_profile(ms_path, flags.nqueries, flags.range_size, flags.from_idx_max, time_usage);
    if(flags.block_size == -1)
        return dispatcher_t::fast_profile(ms_path, flags.nqueries, flags.range_size, flags.from_idx_max, time_usage);
    if(flags.block_size > 0)
        return dispatcher_t::indexed_profile(ms_path, ridx_path, flags.nqueries, flags.range_size, flags.from_idx_max, flags.block_size, time_usage);
    throw string{"bad block_size(" + to_string(flags.block_size) + ") expexting >= 0"};
}


void sdsl_comp(const string& ms_path, const string& ridx_path, rq_dispatcher::counter_t& time_usage, const InputFlags& flags) {
    bool is_rrr = (flags.compression == Compression::none);
    if(flags.block_size == 0){
        if(is_rrr){
            return sdsl_rq_dispatcher::trivial_profile<sdsl::rrr_vector<>>(ms_path, flags.nqueries, flags.range_size, flags.from_idx_max, time_usage);
        } else {
            return sdsl_rq_dispatcher::trivial_profile<sdsl::bit_vector>(ms_path, flags.nqueries, flags.range_size, flags.from_idx_max, time_usage);
        }
    }
    if(flags.block_size > 0){
        if(is_rrr) {
            return sdsl_rq_dispatcher::indexed_profile<sdsl::rrr_vector<>>(ms_path, ridx_path, flags.nqueries, flags.range_size, flags.from_idx_max, flags.block_size, time_usage);
        } else {
            return sdsl_rq_dispatcher::indexed_profile<sdsl::bit_vector>(ms_path, ridx_path, flags.nqueries, flags.range_size, flags.from_idx_max, flags.block_size, time_usage);
        }
    }
    if(flags.block_size == -1){
        if(is_rrr) {
            return sdsl_rq_dispatcher::rrr_fast_profile(ms_path, flags.nqueries, flags.range_size, flags.from_idx_max, flags.compression == Compression::rrr, time_usage);
        } else {
            return sdsl_rq_dispatcher::none_fast_profile(ms_path, flags.nqueries, flags.range_size, flags.from_idx_max, flags.compression == Compression::rrr, time_usage);
        }
    }
    throw string{"bad block_size(" + to_string(flags.block_size) + ") expexting >= 0"};
}

int main(int argc, char **argv) {
    if(argc == 1){
        (cerr << "Run a series of range queries with random end points and report (in csv format) their time in the standard output.\n"
              << "Args:\n"
              << help__ms_path
              << help__ridx_path
              << help__block_size
              << "\t-range_size <positive int>: the size of the range queries\n"
              << "\t-from_max_idx <non-negative int>: range starting point will be greater or equal to this value\n"
              << "\t-niter <positive int>: the number of queries to run\n"
              << "\t-header 1: print a header of the report\n"
              << help__compression
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
    string ms_path = input.getCmdOption("-ms_path");
    string ridx_path = input.getCmdOption("-ridx_path");

    std::srand(123);
    rq_dispatcher::counter_t time_usage;
    switch(flags.compression)
    {
        case Compression::none:
            sdsl_comp(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), time_usage, flags);
            break;
        case Compression::rrr:
            sdsl_comp(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), time_usage, flags);
            break;
        case Compression::rle:
            rle_comp<rle_rq_dispatcher<CSA::RLEVector, CSA::RLEVector::Iterator>>(ms_path, ridx_path, time_usage, flags);
            break;
        case Compression::delta:
            rle_comp<rle_rq_dispatcher<CSA::DeltaVector, CSA::DeltaVector::Iterator>>(ms_path, ridx_path, time_usage, flags);
            break;
        case Compression::nibble:
            rle_comp<rle_rq_dispatcher<CSA::NibbleVector, CSA::NibbleVector::Iterator>>(ms_path, ridx_path, time_usage, flags);
            break;
        case Compression::succint:
            rle_comp<rle_rq_dispatcher<CSA::SuccinctVector, CSA::SuccinctVector::Iterator>>(ms_path, ridx_path, time_usage, flags);
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
