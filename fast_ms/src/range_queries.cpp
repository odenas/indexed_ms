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
#include "fd_ms/p_ms_vector.hpp"
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
    bool check; // can be max/min etc.
    int64_t block_size;
    size_type from_idx, to_idx;
    Compression compression;

    InputFlags() {}

    InputFlags(const InputFlags& f) :
    check{f.check}, from_idx{f.from_idx}, to_idx{f.to_idx},
    compression{f.compression} { }

    InputFlags(bool check, int64_t block_size, size_type from_idx, size_type to_idx) :
    check{check}, block_size{block_size}, from_idx{from_idx}, to_idx{to_idx},
    compression{compression} { }

    InputFlags(OptParser input) :
    check{input.getCmdOption("-check") == "1"}, // check answer
    from_idx{static_cast<size_type> (std::stoi(input.getCmdOption("-from_idx")))},
    to_idx{static_cast<size_type> (std::stoi(input.getCmdOption("-to_idx")))},
    block_size(static_cast<size_type> (std::stoi(input.getCmdOption("-block_size")))) {
        compression = ms_compression::parse_compression(
            input.getCmdOption("-compression")
        );
        if(! valid_range())
            throw string{"Bad range " + range_str()};
    }

    bool valid_range() const {
        return (from_idx < to_idx);
    }

    bool range_out_of_bounds(const size_type max_len) const {
        return (from_idx > max_len) || (to_idx > max_len);
    }

    string range_str() const {
        return string{"[" + to_string(from_idx) + ", " + to_string(to_idx) + ")"};
    }
};

template<typename vec_type, typename it_type>
size_type comp_rle(const string ms_path, const string ridx_path, const InputFlags& flags){
    if(flags.block_size == 0)
        return rle_rq_dispatcher<vec_type, it_type>::trivial(ms_path, flags.from_idx, flags.to_idx, flags.check);
    if(flags.block_size > 0)
        return rle_rq_dispatcher<vec_type, it_type>::indexed(ms_path, ridx_path, flags.from_idx, flags.to_idx, flags.block_size, flags.check);
    return rle_rq_dispatcher<vec_type, it_type>::fast(ms_path, flags.from_idx, flags.to_idx, flags.check);
}


template<typename ms_type, typename ms_sel_1_type>
size_type comp(const string& ms_path, const string& ridx_path, const InputFlags& flags) {
    if(flags.block_size == 0)
        return sdsl_rq_dispatcher<ms_type, ms_sel_1_type>::trivial(ms_path, flags.from_idx, flags.to_idx, flags.check);
    if(flags.block_size > 0)
        return sdsl_rq_dispatcher<ms_type, ms_sel_1_type>::indexed(ms_path, ridx_path, flags.from_idx, flags.to_idx, flags.block_size, flags.check);
    if(flags.block_size == -1)
        return sdsl_rq_dispatcher<ms_type, ms_sel_1_type>::fast(ms_path, flags.from_idx, flags.to_idx, flags.compression == Compression::rrr, flags.check);
    throw string{"bad block_size(" + to_string(flags.block_size) + ") expexting >= 0"};
}

int main(int argc, char **argv) {
    if(argc == 1){
        (cerr << "Answer a range query\n"
              << "Args:\n"
              << help__ms_path
              << help__ridx_path
              << help__from_idx
              << help__to_idx
              << help__block_size
              << help__compression
              << endl);
        exit(0);
    }

    OptParser input(argc, argv);
    string ms_path, ridx_path;

    ms_path = input.getCmdOption("-ms_path");
    ridx_path = input.getCmdOption("-ridx_path");
    try{
        InputFlags flags = InputFlags(input);
        size_type answer = 0;
        if(flags.compression == Compression::none or flags.compression == Compression::rrr){
            answer = comp<sdsl::bit_vector, sdsl::bit_vector::select_1_type>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
        } else {
            switch (flags.compression)
            {
            case Compression::rle:
                answer = comp_rle<CSA::RLEVector, CSA::RLEVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
                break;
            case Compression::delta:
                answer = comp_rle<CSA::DeltaVector, CSA::DeltaVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
                break;
            case Compression::nibble:
                answer = comp_rle<CSA::NibbleVector, CSA::NibbleVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
                break;
            case Compression::succint:
                answer = comp_rle<CSA::SuccinctVector, CSA::SuccinctVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
                break;
            default:
                cerr << "Error." << endl;
                return 1;
            }
        }
        (cout << "[" << flags.from_idx << ", " << flags.to_idx << ")"
          << ": " << answer
          << endl);
    } catch (string s) {
        cerr << s << endl;
        return 1;
    }
}
