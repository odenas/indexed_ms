/*
check that elements of a bit vector interval are equal
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "sdsl/bit_vectors.hpp"
#include "sdsl/util.hpp"

#include "fd_ms/help.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/counter.hpp"
#include "fd_ms/p_ms_vector.hpp"

using namespace fdms;
using namespace std;

typedef StreeOhleb<> cst_t;
typedef typename cst_t::size_type size_type;


class InputFlags {
public:
    size_type start, len;

    InputFlags() { }

    InputFlags(const InputFlags& f) : start{f.start}, len{f.len} { }

    InputFlags(const size_type start, const size_type len) : start{start}, len{len} { }

    InputFlags(OptParser input) :
        start{static_cast<size_type> (std::stoll(input.getCmdOption("-start")))},
        len{static_cast<size_type> (std::stoll(input.getCmdOption("-len")))} {}
};

template<typename ms_type, typename ms_sel_0_type, typename ms_sel_1_type>
int comp(const string ms_path, const InputFlags& flags) {
    ms_type ms;
    cerr << "loading ... ";
    sdsl::load_from_file(ms, ms_path);
    cerr << "DONE" << endl;

    size_type v = ms[flags.start];
    size_type cnt = 0;
    for(size_type j = flags.start; j < flags.start + flags.len; j++){
        if(ms[j] != v){
            cnt += 1;
            v = ms[j];
        }
    }
    cout << flags.start << "," << flags.len << "," << cnt << endl;
}

int main(int argc, char **argv) {
    if(argc == 1){
        (cerr << "Answer a range query\n"
              << "Args:\n"
              << help__ms_path
              << "\t-start <non-negative int>: start of a 0-based half-open interval [from_idx, to_idx)\n"
              << "\t-len <non-negative int>: length of a 0-based half-open interval [from_idx, to_idx)\n"
              << endl);
        exit(0);
    }

    OptParser input(argc, argv);
    string ms_path = input.getCmdOption("-ms_path");

    InputFlags flags;
    try{
        flags = InputFlags(input);
    } catch (string s) {
        cerr << s << endl;
        return 1;
    }

    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);
    return comp<sdsl::bit_vector, sdsl::bit_vector::select_0_type, sdsl::bit_vector::select_1_type>(ms_path, flags);
}