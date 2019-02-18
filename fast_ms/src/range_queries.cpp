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

typedef StreeOhleb<> cst_t;
typedef typename cst_t::node_type node_type;
typedef typename cst_t::size_type size_type;
typedef typename cst_t::char_type char_type;
typedef Counter<size_type> counter_t;

counter_t time_usage{};

class InputFlags {
public:
    bool check; // can be max/min etc.
    size_type block_size;
    size_type from_idx, to_idx;

    InputFlags() {}

    InputFlags(const InputFlags& f) :
    check{f.check}, from_idx{f.from_idx}, to_idx{f.to_idx}{}

    InputFlags(bool check, size_type block_size, size_type from_idx, size_type to_idx) :
    check{check}, block_size{block_size}, from_idx{from_idx}, to_idx{to_idx}
    {}

    InputFlags(OptParser input) :
    check{input.getCmdOption("-check") == "1"}, // check answer
    from_idx{static_cast<size_type> (std::stoi(input.getCmdOption("-from_idx")))},
    to_idx{static_cast<size_type> (std::stoi(input.getCmdOption("-to_idx")))},
    block_size(static_cast<size_type> (std::stoi(input.getCmdOption("-block_size")))) {
    }
};

void comp(const string ms_path, const string ridx_path, const InputFlags& flags) {
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);
    size_type answer = 0;
    if(flags.block_size > 0) {
        sdsl::int_vector<64> ridx;
        sdsl::load_from_file(ridx, ridx_path);
        sdsl::bit_vector::select_1_type ms_sel(&ms);

        partial_sums_vector<size_type> psum(ridx_path, flags.block_size);
        answer = psum.range_sum(ms, ridx, ms_sel, flags.from_idx, flags.to_idx);
    } else {
        sdsl::bit_vector::select_1_type ms_sel(&ms);
        answer = partial_sums_vector<size_type>::trivial_range_sum(ms, ms_sel, flags.from_idx, flags.to_idx);
    }
    if (flags.check) {
        sdsl::bit_vector::select_1_type ms_sel(&ms);
        size_type answer_check = partial_sums_vector<size_type>::trivial_range_sum(ms, ms_sel, flags.from_idx, flags.to_idx);
        if (answer != answer_check) {
            cerr << "answer " << answer << " != expected answer " << answer_check << endl;
            exit(1);
        }
    }
    (cout << "[" << flags.from_idx << ", " << flags.to_idx << ")"
          << " using index " << ridx_path << " of block_size " << flags.block_size
          << ": " << answer
          << endl);
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    string ms_path, ridx_path;

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
        comp("range_queries/m.t.ms", "range_queries/m.t.ms.4.ridx", InputFlags(true, 4, 1, 2));
    }

    ms_path = input.getCmdOption("-ms_path");
    ridx_path = input.getCmdOption("-ridx_path");
    comp(ms_path, ridx_path, InputFlags(input));
    return 0;
}