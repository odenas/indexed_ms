/*
save a block of the given ms vector
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "sdsl/bit_vectors.hpp"
#include "sdsl/util.hpp"

#include "fd_ms/help.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"


using namespace fdms;
using namespace std;

typedef unsigned long long size_type;


class InputFlags {
public:
    size_type start, len;
    bool strict;

    InputFlags() { }

    InputFlags(const InputFlags& f) :
            start{f.start}, len{f.len}, strict{f.strict} { }

    InputFlags(const size_type start, const size_type len, const bool strict) :
            start{start}, len{len}, strict{strict} { }

    InputFlags(OptParser input) :
            strict{input.getCmdOption("-strict") == "1"},
            start{static_cast<size_type> (std::stoll(input.getCmdOption("-from_idx")))},
            len{static_cast<size_type> (std::stoll(input.getCmdOption("-block_size")))} {}
};

void dump_subvector(const string ms_path, const string out_path, const InputFlags& flags){
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);

    if(flags.start >= ms.size())
        throw string{"Out of bounds on ms vector (len " + std::to_string(ms.size()) + "). " +
                     "start=" + std::to_string(flags.start)};
    size_type len = flags.len;
    if(!flags.strict){
        if(flags.start + len > ms.size())
            len = ms.size() - flags.start;
    }

    if(flags.start + len > ms.size())
        throw string{"Out of bounds on vector of length(" +
                     std::to_string(ms.size()) + "). " +
                     "start + len=" + std::to_string(flags.start + len)};

    cerr << "[" << flags.start << ", " << len << ")" << endl;
    sdsl::bit_vector out_ms(len, 0);
    for(size_type i = 0; i < len; i++){
        out_ms[i] = ms[flags.start + i];
    }
    sdsl::store_to_file(out_ms, out_path);
}

int main(int argc, char **argv){
    OptParser input(argc, argv);
    string ms_path;

    if(argc == 1){
        (cerr << "Save a block of the given ms vector\n"
              << "Args:\n"
              << help__ms_path
              << help__from_idx
              << help__block_size
              << "\t-strict <0/1>: if set, exits with error if the block_size is too big\n"
              << "\t-out_path <path to output>: where to write the block\n"
              << endl);
        exit(0);
    }
    InputFlags flags;
    try{
        flags = InputFlags(input);
    }
    catch (string s) {
        cerr << s << endl;
        return 1;
    }
    try{
        dump_subvector(input.getCmdOption("-ms_path"),
                       input.getCmdOption("-out_path"), flags);
    } catch (string s) {
        cerr << "ERROR from dump_subvector():" << endl;
        cerr << s << endl;
    }
}


