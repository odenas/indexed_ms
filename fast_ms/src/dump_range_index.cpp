#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <future>
#include <thread>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>

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
        if (block_size < 0) {
            cerr << "positive block size needed. got: " << block_size << endl;
            exit(1);
        }
    }

public:
    bool time_usage; // can be max/min etc.
    size_t block_size;

    InputFlags() {
    }

    InputFlags(const InputFlags& f) :
    time_usage{f.time_usage}, block_size{f.block_size}
    {
    }

    InputFlags(bool time_, size_t block_size) :
    time_usage{time_}, block_size{block_size}
    {
        check();
    }

    InputFlags(OptParser input) : time_usage{input.getCmdOption("-time_usage") == "1"}
    {
        string bs{input.getCmdOption("-block_size")};
        try{
            block_size = static_cast<size_t>(std::stoi(bs));
        } catch (string m){
            cerr << ("Error: " + m + 
                    "A positive block size is needed (-block_size ). ");
            exit(1);
        }
        check();
    }
};

static void show_MS(const InputSpec ispec, std::ostream& out) {
    buff_vec_t ms(ispec.ms_fname, std::ios::in);
    //sdsl::bit_vector ms;
    //sdsl::load_from_file(ms, ispec.ms_fname);

    size_type k = 0;
    for (size_type i = 0; i < ms.size(); i++) {
        if (ms[i] == 1) {
            out << i - (2 * k) << " ";
            k += 1;
        }
    }
}

void comp(const string ms_path, counter_t& time_usage, InputFlags& flags) {
    auto comp_start = timer::now();
    buff_vec_t ms(ms_path, std::ios::in);
    sdsl::int_vector_buffer<64> out_vec(
        ms_path + "." + to_string(flags.block_size) + ".range",
        std::ios::out
    );

    size_type one_cnt = 0, out_idx = 0, ms_value = 0;
    for(size_type ms_idx = 0; ms_idx < ms.size(); ms_idx++){
        if(ms[ms_idx] == 1){
            ms_value += ms_idx - (2 * one_cnt);
            one_cnt += 1;
        }
        if(ms_idx % flags.block_size == 0){
            out_vec[out_idx++] = ms_value;
            cout << "out_vec[" << out_idx << "] = " << out_vec[out_idx - 1] << endl;
        }
    }
    time_usage.register_now("comp_total", comp_start);
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    string ms_path;
    InputFlags flags;
    counter_t time_usage{};

    if (argc == 1) {
        ms_path = {"/home/brt/code/matching_statistics/indexed_ms/fast_ms/tests/"};
        flags = InputFlags(
                false, // time
                4 // block_size
                );
    } else {
        ms_path = input.getCmdOption("-ms_path");
        flags = InputFlags(input);
    }

    cout << "ms_path: " << ms_path << endl;

    auto comp_start = timer::now();
    comp(ms_path, time_usage, flags);
    time_usage.register_now("comp_total", comp_start);

    if (flags.time_usage) {
        cerr << "dumping reports ..." << endl;
        for (auto item : time_usage.reg)
            (cout << item.first << "," << item.second << endl);
    }
    return 0;
}
