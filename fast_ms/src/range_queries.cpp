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


class InputFlags {
private:
    void check() const {
        if(from >= to){
            cerr << "empty interval [" << from << ", " << to << ")" << endl;
            exit(1);
        }
    }

public:
    bool time_usage, avg; // can be max/min etc.
    size_type from, to, block_size;

    InputFlags() { }

    InputFlags(const InputFlags& f) : 
    time_usage{f.time_usage}, avg{f.avg}, from{f.from}, to{f.to} { }

    InputFlags(bool time_, bool avg, size_t from, size_t to) : 
    time_usage{time_}, avg{avg}, from{from}, to{to}
    {
        check();
    }

    InputFlags(OptParser input) :
    time_usage{input.getCmdOption("-time_usage") == "1"}, // time usage
    avg{input.getCmdOption("-avg") == "1"},               // average matching statistics
    from{static_cast<size_type> (std::stoi(input.getCmdOption("-from")))},
    to{static_cast<size_type> (std::stoi(input.getCmdOption("-to")))},
    block_size(static_cast<size_type>(std::stoi(input.getCmdOption("-block_size")))){
        check(); 
    }
};

size_type sum_ms_prefix(const sdsl::bit_vector& ms, const size_type to_ms_idx){
    size_type prev_ms = 1, sum_ms = 0;
    size_type cnt1 = 0, cnt0 = 0;
    size_type i = 0;
    while(i <= to_ms_idx){
        if(ms[i++] == 1){
            //(cout << "MS[" << cnt1 - 1 << "] = " << prev_ms << ", SUM = " << sum_ms << endl);
            prev_ms += (cnt0 - 1);
            sum_ms += prev_ms;
            cnt1 += 1;
            cnt0 = 0;
        }
        else{
            cnt0 += 1;
        }
    }
    //(cout << "MS[" << cnt1 - 1 << "] = " << prev_ms << ", SUM = " << sum_ms << endl);
    return sum_ms;
}

void comp(const string ms_path, const string ridx_path, InputFlags& flags) {
    auto comp_start = timer::now();
    /*
    sdsl::bit_vector ms = {
        0, 0, 0, 1,
        0, 0, 1,
        1,
        1,
        0, 0, 1,
        1,
        0, 0, 1,
        0, 1,
        0, 0, 1,
        1,
        1,
        1,
    };
    sdsl::int_vector<64> ridx = {3, 10, 15, 20, 23, 33};
    */

    sdsl::bit_vector ms; sdsl::load_from_file(ms, ms_path);
    //cerr << "loaded ms bit-vector (size " << ms.size() << ") from " << ms_path << endl;
    //cerr << "entries:" << endl;
    //for(int i=0; i < ms.size(); i++)
    //   cerr << i << " ) " << ms[i] << endl;
    
    sdsl::int_vector_buffer<64> ridx(
        ms_path + "." + to_string(flags.block_size) + ".range",
        std::ios::in
    );
    //(cerr << "loaded msidx bit-vector (size " << ridx.size() << ") from " << 
    //        ms_path + "." + to_string(flags.block_size) + ".range" << endl);
    //cerr << "entries:" << endl;
    //for(int i=0; i < ridx.size(); i++)
    //    cerr << i << " ) " << ridx[i] << endl;

    sdsl::bit_vector::select_1_type ms_sel(&ms);
    sdsl::bit_vector::rank_1_type ms_rank(&ms);

    // index of last term of sum
    size_type int_ms_idx = flags.to - 1;
    size_type bit_ms_idx = ms_sel(int_ms_idx + 1);
    size_type block_idx = bit_ms_idx / flags.block_size;
    
    size_type sum_ms = 0;  // to be subtracted from ridx[block_idx]
    {
        size_type prev_ms = bit_ms_idx - (2 * int_ms_idx); // needed for 1st term beyond the sum
        size_type end_block_nones = ms_rank((block_idx + 1) * flags.block_size);    
        size_type nzeros = 0, nones = int_ms_idx;

        // loop from bit_ms_idx + 1 to the end of the block
        for(size_type i = bit_ms_idx + 1; i < (block_idx + 1) * flags.block_size; i++){
            if(nones >= end_block_nones)
                break;

            if(ms[i] == 1){
                size_type cur_ms = (prev_ms + nzeros - 1);
                sum_ms += cur_ms;  // since MS_i - MS_{i-1} + 1 = nzeros
                nones += 1;
                prev_ms = cur_ms;
            } else {
                nzeros += 1;
            }
        }
    }
    (cerr 
          << "bit_ms_idx = " << bit_ms_idx << " "
          << "block_idx = " << block_idx << " "
          << "sum_ms = " << sum_ms << " "
          << endl);
    cout << ridx[block_idx] - sum_ms << endl;
}


int main(int argc, char **argv) {
    OptParser input(argc, argv);
    string ms_path, ridx_path;
    InputFlags flags;
    counter_t time_usage{};

    if (argc == 1) {
        exit(1);
        ms_path = {"/home/brt/code/matching_statistics/indexed_ms/fast_ms/tests/"};
        flags = InputFlags(
                false, // time
                false, // avg
                0, 3   // range
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
