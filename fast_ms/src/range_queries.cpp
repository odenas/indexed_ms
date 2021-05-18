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
#include "fd_ms/partial_max_vector.hpp"
#include "fd_ms/range_query.hpp"
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
typedef rq_result<size_type> rqres_t;

class InputFlags {
    static RangeOperation parse_operation(const string c_str){
        std::map<RangeOperation, string> a2s = {
            {RangeOperation::r_sum, "sum"},
            {RangeOperation::r_max, "max"}
        };

        if(c_str == "0" or c_str == "sum")
            return RangeOperation::r_sum;
        for(auto item: a2s){
            if(item.second == c_str)
                return item.first;
        }
        throw string{"bad operation string: " + c_str};
    }

    static RangeAlgorithm parse_algo(const string c_str){
        std::map<RangeAlgorithm, string> a2s = {
            {RangeAlgorithm::trivial, ".t"},
            {RangeAlgorithm::djamal, ".d"}
        };

        if(c_str == "0" or c_str == "none")
            return RangeAlgorithm::none;
        for(auto item: a2s){
            if(item.second == ("." + c_str))
                return item.first;
        }
        throw string{"bad compression string: " + c_str};
    }

public:
    int64_t block_size;
    size_type from_idx, to_idx;
    Compression compression;
    RangeAlgorithm algo;
    RangeOperation op;

    InputFlags() {}

    InputFlags(const InputFlags& f) :
    from_idx{f.from_idx}, to_idx{f.to_idx},
    compression{f.compression},
        algo{f.algo}, op{f.op} { }

    InputFlags(OptParser input) :
    from_idx{static_cast<size_type> (std::stoi(input.getCmdOption("-from_idx")))},
    to_idx{static_cast<size_type> (std::stoi(input.getCmdOption("-to_idx")))},
    block_size(static_cast<size_type> (std::stoi(input.getCmdOption("-block_size")))) {
        if(block_size < 0)
            throw string{"Bad block size."};

        compression = ms_compression::parse_compression(
            input.getCmdOption("-compression")
        );
        if(! valid_range())
            throw string{"Bad range " + range_str()};

        algo = parse_algo(input.getCmdOption("-algo"));
        if(algo == RangeAlgorithm::none)
            throw string{"Expecting an RangeAlgorithm."};
        op = parse_operation(input.getCmdOption("-op"));
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
rqres_t comp_rle(const string ms_path, const string ridx_path, const InputFlags& flags){
    if(flags.block_size == 0)
        return rle_rq_dispatcher<vec_type, it_type>::noindex(ms_path, flags.from_idx, flags.to_idx,
            flags.algo, flags.op);
    return rle_rq_dispatcher<vec_type, it_type>::indexed(ms_path, ridx_path,
        flags.from_idx, flags.to_idx, flags.block_size, flags.algo, flags.op);
}


rqres_t comp_sdsl(const string& ms_path, const string& ridx_path, const InputFlags& flags) {
    bool is_rrr = (flags.compression == Compression::rrr);
    if(flags.block_size == 0){
        if(is_rrr)
            return sdsl_rq_dispatcher::rrr_noindex(
                ms_path, flags.from_idx, flags.to_idx, flags.algo, flags.op);
        return sdsl_rq_dispatcher::none_noindex(
            ms_path, flags.from_idx, flags.to_idx, flags.algo, flags.op);
    }
    if(flags.block_size > 0){
        if(is_rrr) {
            return sdsl_rq_dispatcher::rrr_indexed(ms_path, ridx_path,
                flags.from_idx, flags.to_idx, flags.block_size, flags.algo, flags.op);
        }
        return sdsl_rq_dispatcher::none_indexed(ms_path, ridx_path,
            flags.from_idx, flags.to_idx, flags.block_size, flags.algo, flags.op);
    }
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
              << help__algo
              << help__rangeop
              << endl);
        exit(0);
    }

    OptParser input(argc, argv);
    try{
        InputFlags flags = InputFlags(input);
        rqres_t answer;
        switch (flags.compression)
        {
            case Compression::none:
                answer = comp_sdsl(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
                break;
            case Compression::rrr:
                answer = comp_sdsl(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
                break;
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
        (cout << "[" << flags.from_idx << ", " << flags.to_idx << ")"
          << ": " << answer.value
          << endl);
    } catch (string s) {
        cerr << s << endl;
        return 1;
    }
}
