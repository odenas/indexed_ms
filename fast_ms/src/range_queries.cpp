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
    size_t from, to;

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
    avg{input.getCmdOption("-avg") == "1"},                // average matching statistics
    from{static_cast<size_t> (std::stoi(input.getCmdOption("-from")))},
    to{static_cast<size_t> (std::stoi(input.getCmdOption("-to")))}
    { 
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

void comp(const string ms_path, counter_t& time_usage, InputFlags& flags) {
    auto comp_start = timer::now();
    sdsl::bit_vector ms = {0, 0, 1,    //2
                           1,          //1
                           0, 0, 0, 1, //3
                           0, 1,       //3
                           1,          //2
                           1};         //1
    
    
    //sdsl::bit_vector ms; sdsl::load_from_file(ms, ms_path);
    //cerr << "loaded ms bit-vector (size " << ms.size() << ") from " << ms_path << endl;

    sdsl::bit_vector::select_1_type ms_sel(&ms);

    size_type sum_ms = sum_ms_prefix(ms, ms_sel(flags.to));
    if(flags.from)
        sum_ms -= sum_ms_prefix(ms, ms_sel(flags.from));

    time_usage.register_now("comp_total", comp_start);

    (cout << "AVG(MS[" << flags.from << ", " << flags.to - 1<< "]) = " 
            << sum_ms << " / " << (flags.to - flags.from) << endl);
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
                false, // avg
                0, 3   // range
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
