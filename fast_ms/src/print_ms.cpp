/*
print sections of an MS bit vector
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

using namespace fdms;
using namespace std;

typedef unsigned long long size_type;


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
    sdsl::load_from_file(ms, ms_path);

    for(size_type j = flags.start; j < flags.start + flags.len; j++){
        cout << static_cast<int>(ms[j]);
        if(j >= ms.size()){
            cerr << "reached the end at " << j << endl;
            break;
        }
    }
    cout << endl;
}

int main(int argc, char **argv) {
    if(argc == 1){
        (cerr << "Print a section of the ms bit-vector\n"
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
    return comp<sdsl::bit_vector, sdsl::bit_vector::select_0_type, sdsl::bit_vector::select_1_type>(ms_path, flags);
}
