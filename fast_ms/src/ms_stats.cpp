/*
print stats of an MS bit vector
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

using namespace fdms;
using namespace std;

typedef unsigned long long size_type;
typedef sdsl::int_vector_buffer<1> buff_vec_t;


class InputFlags {
public:
    size_type start, len;

    InputFlags() { }

    InputFlags(const InputFlags& f) :
        start{f.start}, len{f.len}
    { }

    InputFlags(const size_type start, const size_type len, const bool int_format) :
        start{start}, len{len}
    { }

    InputFlags(OptParser input) :
        start{static_cast<size_type> (std::stoll(input.getCmdOption("-start")))},
        len{static_cast<size_type> (std::stoll(input.getCmdOption("-len")))} {}
};

int comp(const string ms_path, const InputFlags& flags) {
    buff_vec_t ms(ms_path, std::ios::in);

    size_type end = flags.start + flags.len;
    if(flags.len == 0)
        end = ms.size();

    size_type sum = 0;

    for(size_type j = flags.start; j < end; j++){
        sum += ms[j];
        if(j >= ms.size()){
            cerr << "reached the end at " << j << endl;
            break;
        }
    }
    cout << "filename:    " << ms_path << endl;
    cout << "total size:  " << ms.size() << endl;
    cout << "range:       [" << flags.start << ", " << end << ")" << endl;
    cout << "range size:  " << end - flags.start << endl;
    cout << "0s in range: " << end - flags.start - sum << endl;
    cout << "1s in range: " << sum << endl;
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

