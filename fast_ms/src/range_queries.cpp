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

using namespace std;
using namespace fdms;

typedef StreeOhleb<> cst_t;
typedef typename cst_t::node_type node_type;
typedef typename cst_t::size_type size_type;
typedef typename cst_t::char_type char_type;
typedef Counter<size_type> counter_t;
typedef sdsl::int_vector_buffer<1> buff_vec_t;
std::vector<size_type> end_points = {
    1043086, 986408, 1007099, 1186359, 
    1130095, 1093937, 957132, 1124650, 
    1111185, 848593, 1168285, 1110883, 
    1083327, 803228, 1061229, 1070160, 
    939141, 977511, 1108004, 1185944};
#define LOOPS 5

counter_t time_usage{};


class InputFlags {
public:
    bool time_usage, check; // can be max/min etc.
    size_type from, to, block_size;

    InputFlags() {
    }

    InputFlags(const InputFlags& f) :
    time_usage{f.time_usage}, check{f.check}, from{f.from}, to{f.to}
    {
    }

    InputFlags(bool time_, bool check, size_t from, size_t to) :
    time_usage{time_}, check{check}, from{from}, to{to}
    {
        if (from >= to) {
            cerr << "empty interval [" << from << ", " << to << ")" << endl;
            exit(1);
        }
    }

    InputFlags(OptParser input) :
    time_usage{input.getCmdOption("-time_usage") == "1"}, // time usage
    check{input.getCmdOption("-check") == "1"}, // check answer
    from{static_cast<size_type> (std::stoi(input.getCmdOption("-from")))},
    to{static_cast<size_type> (std::stoi(input.getCmdOption("-to")))}
,
    block_size(static_cast<size_type> (std::stoi(input.getCmdOption("-block_size")))) {
    }
};

size_type sum_ms_prefix(const sdsl::bit_vector& ms, const size_type to_ms_idx) {
    //cerr << "[trivial] " << to_ms_idx << endl;
    
    size_type prev_ms = 1, cur_ms = 0, sum_ms = 0;
    size_type cnt1 = 0, cnt0 = 0, i = 0;
    while (cnt1 < to_ms_idx) {
        if (ms[i] == 1) {
            //(cout << "MS[" << cnt1 - 1 << "] = " << prev_ms << ", SUM = " << sum_ms << endl);
            cur_ms = prev_ms + cnt0 - 1;
            sum_ms += cur_ms;
            prev_ms = cur_ms;
            cnt0 = 0;
            cnt1 += 1;
        } else {
            cnt0 += 1;
        }
        i += 1;
    }
    return sum_ms;
}

size_type sum_ms_prefix(sdsl::bit_vector& ms, sdsl::int_vector<64>& ridx,  
        sdsl::bit_vector::select_1_type ms_sel, sdsl::bit_vector::rank_1_type ms_rank,
        const size_type to_ms_idx, const InputFlags& flags) {
    //cerr << "[indexed (" << flags.block_size << ")] " << to_ms_idx << endl;

    // index of last term of sum
    size_type int_ms_idx = to_ms_idx;
    size_type bit_ms_idx = ms_sel(int_ms_idx + 1);
    size_type block_idx = bit_ms_idx / flags.block_size;

    size_type sum_ms = 0; // to be subtracted from ridx[block_idx]
    {
        size_type prev_ms = bit_ms_idx - (2 * int_ms_idx); // needed for 1st term beyond the sum
        size_type end_block_nones = ms_rank((block_idx + 1) * flags.block_size);
        size_type nzeros = 0, nones = int_ms_idx;

        // loop from bit_ms_idx + 1 to the end of the block
        for (size_type i = bit_ms_idx + 1; i < (block_idx + 1) * flags.block_size; i++) {
            if (nones >= end_block_nones)
                break;

            if (ms[i] == 1) {
                size_type cur_ms = (prev_ms + nzeros - 1);
                sum_ms += cur_ms; // since MS_i - MS_{i-1} + 1 = nzeros
                nones += 1;
                prev_ms = cur_ms;
            } else {
                nzeros += 1;
            }
        }
    }
    size_type answer = ridx[block_idx] - sum_ms;
    return answer;
}

void comp(const string ms_path, const string ridx_path, const InputFlags& flags) {
    auto comp_start = timer::now();
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);
    time_usage.register_now("load_ms", comp_start);

    size_type answer = 0;
    if(flags.block_size == 0){
        time_usage.register_now("load_partial_sums", timer::now());
        time_usage.register_now("rank_sel_init", timer::now());

        auto comp_start = timer::now();
        for (int k = 0; k < end_points.size() * LOOPS; k++){
            answer = sum_ms_prefix(ms, flags.to + end_points[k % end_points.size()]);
        }
        time_usage.register_now("algorithm", comp_start);

    } else {
        auto ds_start = timer::now(); // should be 0 since it is just the buffer
        //sdsl::int_vector_buffer<64> ridx(ridx_path, std::ios::in);
        sdsl::int_vector<64> ridx;
        sdsl::load_from_file(ridx, ridx_path);
        time_usage.register_now("load_partial_sums", ds_start);

        ds_start = timer::now();
        sdsl::bit_vector::select_1_type ms_sel(&ms);
        sdsl::bit_vector::rank_1_type ms_rank(&ms);
        time_usage.register_now("rank_sel_init", ds_start);

        auto comp_start = timer::now();
        for (int k = 0; k < end_points.size() * LOOPS; k++){
            answer = sum_ms_prefix(ms, ridx, 
                    ms_sel, ms_rank,
                    flags.to + end_points[k % end_points.size()], flags);
        }
        time_usage.register_now("algorithm", comp_start);
    }

    comp_start = timer::now();
    if (flags.check) {
        size_type answer_check = sum_ms_prefix(ms, flags.to);
        if (answer != answer_check) {
            cerr << "answer " << answer << " != expected answer " << answer_check << endl;
            exit(1);
        }
    }
    time_usage.register_now("check", comp_start);
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    string ms_path, ridx_path;
    InputFlags flags;

    if (argc == 1) {
        exit(1);
        ms_path = {"/home/brt/code/matching_statistics/indexed_ms/fast_ms/tests/"};
        flags = InputFlags(
                false, // time
                false, // avg
                0, 3 // range
                );
    } else {
        ms_path = input.getCmdOption("-ms_path");
        ridx_path = input.getCmdOption("-ridx_path");
        flags = InputFlags(input);
    }
    auto comp_start = timer::now();
    comp(ms_path, ridx_path, flags);
    time_usage.register_now("comp_total", comp_start);



    if (flags.time_usage) {
        cerr << "dumping reports ..." << endl;
        for (auto item : time_usage.reg)
            (cout << item.first << "," << item.second << endl);
    }
    
    return 0;
}
