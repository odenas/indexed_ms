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
#include "fd_ms/slices.hpp"


using namespace fdms;
using namespace std;

typedef unsigned long long size_type;


class InputFlags {
public:
    size_type len;

    InputFlags() { }

    InputFlags(const InputFlags& f) : len{f.len} { }

    InputFlags(const size_type len) : len{len} { }

    InputFlags(OptParser input) :
            len{static_cast<size_type> (std::stoll(input.getCmdOption("-block_size")))} {}
};

void dump_subvector(const sdsl::bit_vector& ms, const size_type start, const size_type end, const string out_path){
    if(start >= ms.size())
        throw string{"Out of bounds on ms vector (len " + std::to_string(ms.size()) + "). " +
                     "start=" + std::to_string(start)};

    size_type w = 64, len = end - start;
    cerr << "[" << start << ", " << end << ")" << endl;
    sdsl::bit_vector out_ms(end - start, 0);

    size_type i = 0;
    while(i + w < len){
        out_ms.set_int(i, ms.get_int(start + i), w);
        i += w;
    }
    while(i < len){
        out_ms[i] = ms[start + i];
        i += 1;
    }

    // check
    for(size_type i = start; i < end; i++){
        if(out_ms[i - start] != ms[i])
            throw string{"Check error"};
    }
    //for(size_type i = start; i < end; i++){
    //    out_ms[i - start] = ms[i];
    //}
    sdsl::store_to_file(out_ms, out_path);
}

int main(int argc, char **argv){
    OptParser input(argc, argv);
    string ms_path;

    if(argc == 1){
        (cerr << "Chop the the given ms vector in blocks of given size.\n"
              << "Args:\n"
              << help__ms_path
              << help__block_size
              << "\t-out_prefix <prefix of outputs>: files will be <prefix>_<start>_<block_size>_<suffix>\n"
              << "\t-out_suffix <suffix of outputs>: files will be <prefix>_<start>_<block_size>_<suffix>\n"
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
        sdsl::bit_vector ms;
        sdsl::load_from_file(ms, input.getCmdOption("-ms_path"));
        string prefix = input.getCmdOption("-out_prefix");
        string suffix = input.getCmdOption("-out_suffix");

        BlockSlices<size_type> slices(ms.size(), flags.len);
        for(int slice_idx = 0; slice_idx < slices.slices.size(); slice_idx++){
            string out_path = (prefix +
                               "_" + to_string(slices[slice_idx].first) +
                               "_" + to_string(slices[slice_idx].second) +
                               suffix);
            dump_subvector(ms, slices[slice_idx].first, slices[slice_idx].second, out_path);
        }
    } catch (string s) {
        cerr << "ERROR from dump_subvector():" << endl;
        cerr << s << endl;
    }
}


