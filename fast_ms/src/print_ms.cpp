/*
print sections of an MS bit vector
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "sdsl/vectors.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/util.hpp"

#include "fd_ms/help.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/counter.hpp"

using namespace fdms;
using namespace std;

typedef unsigned long long size_type;
typedef sdsl::int_vector_buffer<1> buff_vec_t;


class InputFlags {
public:
    size_type start, len;
    bool int_format;

    InputFlags() { }

    InputFlags(const InputFlags& f) :
        start{f.start}, len{f.len}, int_format{f.int_format}
    { }

    InputFlags(const size_type start, const size_type len, const bool int_format) :
        start{start}, len{len}, int_format{int_format}
    { }

    InputFlags(OptParser input) :
        int_format{input.getCmdOption("-int_format") == "1"},
        start{static_cast<size_type> (std::stoll(input.getCmdOption("-start")))},
        len{static_cast<size_type> (std::stoll(input.getCmdOption("-len")))} {}
};

int comp(const string ms_path, const InputFlags& flags) {
    buff_vec_t ms(ms_path, std::ios::in);
    //sdsl::bit_vector ms;
    //sdsl::load_from_file(ms, ms_path);

    for(size_type j = flags.start; j < flags.start + flags.len; j++){
        cout << static_cast<int>(ms[j]);
        if(j >= ms.size()){
            cerr << "reached the end at " << j << endl;
            break;
        }
    }
    cout << endl;
    return 0;
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
    return comp(ms_path, flags);
}
