/*
split an ms vector into blocks of given type
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "sdsl/bit_vectors.hpp"
#include "sdsl/util.hpp"

#include "fd_ms/help.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/ms_block_types.hpp"


using namespace fdms;
using namespace std;

typedef typename ms_blocks::size_type size_type;


class InputFlags {
public:
    string block_t, prefix, suffix;
    size_type len;
    bool check;

    InputFlags() { }

    InputFlags(const InputFlags& f) :
        block_t{f.block_t},
        prefix{f.prefix}, suffix{f.suffix},
        len{f.len}, check{f.check} { }

    InputFlags(OptParser input) :
            len{static_cast<size_type> (std::stoll(input.getCmdOption("-len")))},
            check{input.getCmdOption("-no_check") == "0"}
    {
        block_t = input.getCmdOption("-block_type");
        prefix = input.getCmdOption("-out_prefix");
        suffix = input.getCmdOption("-out_suffix");
    }
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

    sdsl::store_to_file(out_ms, out_path);
}

template<typename block_type_cls>
int dump_generic(const sdsl::bit_vector& ms, const InputFlags flags){
    block_type_cls blocks(ms, flags.len);

    for(size_type slice_idx = 0; slice_idx < blocks.slices.size(); slice_idx++){
        string out_path = blocks.make_out_fname(flags.prefix, flags.suffix, slice_idx);
        dump_subvector(ms, blocks.slices[slice_idx].first, blocks.slices[slice_idx].second, out_path);

        if(flags.check){
            try{ blocks.check(out_path, slice_idx); }
            catch (string s) { throw s; }
        }
    }

    return 0;
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    string ms_path;

    if(argc == 1){
        (cerr << "Chop the the given ms vector in blocks of given size or type.\n"
              << "Args:\n"
              << help__ms_path
              << "\t-len <positive int>: Length of blocks \n"
              << "\t-block_type <block type>: One of: \n"
              << "\t\tfixed_size - blocks of size <block_size>; \n"
              << "\t\tconstant_ones - blocks with <block_size> ones; \n"
              << "\t\tzeros_and_ones - of the type 0^i1^j of size at least <block_size>.\n"
              << "\t-out_prefix <prefix of outputs>: files will be <prefix>_<start>_<block_size><suffix>\n"
              << "\t-out_suffix <suffix of outputs>: files will be <prefix>_<start>_<block_size><suffix>\n"
              << "\t-no_check <0/1>: disable correctness check for each blocks\n"
              << endl);
        exit(0);
    }
    InputFlags flags = InputFlags(input);
    try{
        sdsl::bit_vector ms;
        sdsl::load_from_file(ms, input.getCmdOption("-ms_path"));

        if(flags.block_t == "fixed_size")
            return dump_generic<ms_blocks>(ms, flags);
        else if(flags.block_t == "constant_ones")
            return dump_generic<const_ones_ms_blocks>(ms, flags);
        else if(flags.block_t == "zeros_and_ones")
            return dump_generic<zo_patt_ms_blocks>(ms, flags);
        else
            throw string("Unknown block type" + flags.block_t);
    } catch (string s) {
        cerr << "ERROR:" << endl;
        cerr << s << endl;
    }
}


