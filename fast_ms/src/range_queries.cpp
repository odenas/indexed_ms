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
    bool check; // can be max/min etc.
    size_type block_size;
    size_type from_idx, to_idx;
    Compression compression;

    InputFlags() {}

    InputFlags(const InputFlags& f) :
    check{f.check}, from_idx{f.from_idx}, to_idx{f.to_idx},
    compression{f.compression} { }

    InputFlags(bool check, size_type block_size, size_type from_idx, size_type to_idx) :
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
    }
};

template<typename vec_type, typename it_type>
int comp1(const string ms_path, const string ridx_path, const InputFlags& flags){
    size_type answer = 0, answer_check = 0;

    std::ifstream in{ms_path, std::ios::binary};
    vec_type ms(in);
    it_type* it = new it_type(ms);

    if(flags.block_size == 0 or flags.check)
        answer_check = partial_sums_vector1<vec_type, it_type, size_type>::trivial_range_sum(ms, it, flags.from_idx, flags.to_idx);

    if(flags.block_size > 0){
        sdsl::int_vector<64> ridx;
        sdsl::load_from_file(ridx, ridx_path);

        partial_sums_vector1<vec_type, it_type, size_type> psum(ridx_path, flags.block_size);
        answer = psum.range_sum(ms, ridx, it, flags.from_idx, flags.to_idx);
        if(flags.check)
            answer_check = partial_sums_vector1<vec_type, it_type, size_type>::trivial_range_sum(ms, it, flags.from_idx, flags.to_idx);
    } else {
        answer = partial_sums_vector1<vec_type, it_type, size_type>::trivial_range_sum(ms, it, flags.from_idx, flags.to_idx);
        if(flags.check)
            answer_check = answer;
    }
    if (flags.check && answer != answer_check) {
        cerr << "answer " << answer << " != expected answer " << answer_check << endl;
        exit(1);
    }
    (cout << "[" << flags.from_idx << ", " << flags.to_idx << ")"
          << " using index " << ridx_path << " of block_size " << flags.block_size
          << ": " << answer
          << endl);
    return 0;
}

template<typename ms_type, typename ms_sel_1_type>
int comp(const string ms_path, const string ridx_path, const InputFlags& flags) {
    ms_type ms;
    sdsl::load_from_file(ms, ms_path);
    size_type answer = 0, answer_check = 0;
    ms_sel_1_type ms_sel(&ms);

    if(flags.block_size > 0) {
        sdsl::int_vector<64> ridx;
        sdsl::load_from_file(ridx, ridx_path);

        partial_sums_vector<size_type, ms_type, ms_sel_1_type> psum(ridx_path, flags.block_size);
        answer = psum.range_sum(ms, ridx, ms_sel, flags.from_idx, flags.to_idx);
    } else {
        answer = partial_sums_vector<size_type, ms_type, ms_sel_1_type>::trivial_range_sum(ms, ms_sel, flags.from_idx, flags.to_idx);
    }
    if (flags.check) {
        answer_check = partial_sums_vector<size_type, ms_type, ms_sel_1_type>::trivial_range_sum(ms, ms_sel, flags.from_idx, flags.to_idx);
        if (answer != answer_check) {
            cerr << "answer " << answer << " != expected answer " << answer_check << endl;
            exit(1);
        }
    }
    (cout << "[" << flags.from_idx << ", " << flags.to_idx << ")"
          << " using index " << ridx_path << " of block_size " << flags.block_size
          << ": " << answer
          << endl);
    return 0;
}

int main(int argc, char **argv) {
    if(argc == 1){
        (cerr << "Answer a range query\n"
              << "Args:\n"
              << help__ms_path
              << "\t-ridx_path <path to binary file>: the range query index; only used if block_size > 0\n"
              << "\t-from_idx <non-negative int>: start of a 0-based half-open interval [from_idx, to_idx)\n"
              << "\t-to_idx <non-negative int>: end of a 0-based half-open interval [from_idx, to_idx)\n"
              << "\t-block_size <non-negative int>: range query index block size. If 0, do not use index.\n"
              << endl);
        exit(0);
    }

    OptParser input(argc, argv);
    string ms_path, ridx_path;
    InputFlags flags;

    ms_path = input.getCmdOption("-ms_path");
    ridx_path = input.getCmdOption("-ridx_path");
    try{
        flags = InputFlags(input);
    } catch (string s) {
        cerr << s << endl;
        return 1;
    }

    switch(flags.compression)
    {
        case Compression::none:
            return comp<sdsl::bit_vector, sdsl::bit_vector::select_1_type>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
        case Compression::hyb:
            return comp<sdsl::hyb_vector<>, sdsl::hyb_vector<>::select_1_type>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
        case Compression::rrr:
            return comp<sdsl::rrr_vector<>, sdsl::rrr_vector<>::select_1_type>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
        case Compression::rle:
            return comp1<CSA::RLEVector, CSA::RLEVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
        case Compression::delta:
            return comp1<CSA::DeltaVector, CSA::DeltaVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
        case Compression::nibble:
            return comp1<CSA::NibbleVector, CSA::NibbleVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
        case Compression::succint:
            return comp1<CSA::SuccinctVector, CSA::SuccinctVector::Iterator>(input.getCmdOption("-ms_path"), input.getCmdOption("-ridx_path"), flags);
        default:
            cerr << "Error." << endl;
            return 1;
    }
}