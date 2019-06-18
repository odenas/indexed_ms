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

#include "rlcsa/bits/bitvector.h"
#include "rlcsa/bits/rlevector.h"
#include "rlcsa/bits/deltavector.h"
#include "rlcsa/bits/succinctvector.h"
#include "rlcsa/bits/nibblevector.h"

using namespace fdms;
using namespace std;

typedef typename ms_blocks::size_type size_type;
typedef typename ms_blocks::pair_t pair_t;


class InputFlags {
public:
    string block_t, prefix, suffix, compr;
    size_type len;
    bool check;

    InputFlags() { }

    InputFlags(const InputFlags& f) :
        block_t{f.block_t}, compr{f.compr},
        prefix{f.prefix}, suffix{f.suffix},
        len{f.len}, check{f.check} { }

    InputFlags(OptParser input) :
            len{static_cast<size_type> (std::stoll(input.getCmdOption("-len")))},
            check{input.getCmdOption("-no_check") == "0"}
    {
        block_t = input.getCmdOption("-block_type");
        prefix = input.getCmdOption("-out_prefix");
        suffix = input.getCmdOption("-out_suffix");
        compr = input.getCmdOption("-compr");
        cerr << block_t << endl;
    }
};

template<typename enc_type>
size_type fill_encoder(sdsl::bit_vector ms, enc_type& encoder){
    size_type no = 0, i = 0, n_runs = 0;
    while(i < ms.size()){
        if(ms[i] == 1)
            no += 1;
        else{
            assert (i >= no);
            if (no > 0){
                n_runs += 1;
                encoder.addRun(i - no, no);
                no = 0;
            }
        }
        i += 1;
    }
    if(no > 0){
        encoder.addRun(i - no, no);
        n_runs += 1;
    }
    encoder.flush();
    return n_runs;
}

template<typename vec_type, typename enc_type>
int comp1(const sdsl::bit_vector& ms){
    enc_type encoder(32);
    size_type n_runs = fill_encoder<enc_type>(ms, encoder);
    vec_type c_ms(encoder, ms.size());
    return 0;
}

void dump_subvector(const sdsl::bit_vector& ms,
        const pair_t range, const string out_path, const InputFlags& flags){

    size_type start = range.first;
    size_type end = range.second;

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
    //sdsl::store_to_file(out_ms, out_path);

    if(flags.compr == "none"){
        sdsl::bit_vector::select_1_type a(&out_ms);
    }
    else if (flags.compr == "baseline"){
        ;
    }
    else if(flags.compr == "rrr"){
        sdsl::rrr_vector<>(out_ms);
        sdsl::rrr_vector<>::select_1_type a(&out_ms);
    }
    // RLCSA
    else if (flags.compr == "rle"){
        comp1<CSA::RLEVector, CSA::RLEEncoder>(out_ms);
    }
    else
        throw string{"Bad compression string: " + flags.compr};
}

template<typename block_type_cls>
int dump_generic(const string ms_path, const InputFlags& flags){
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);

    block_type_cls blocks(ms, flags.len);

    for(size_type slice_idx = 0; slice_idx < blocks.slices.size(); slice_idx++){
        string out_path = blocks.make_out_fname(flags.prefix, flags.suffix, slice_idx);
        dump_subvector(ms, blocks.slices[slice_idx], out_path, flags);

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
              << "\t-compr <compression type>: One of: rrr, rle, none, or baseline.\n"
              << "\t-out_prefix <prefix of outputs>: files will be <prefix>_<start>_<block_size><suffix>\n"
              << "\t-out_suffix <suffix of outputs>: files will be <prefix>_<start>_<block_size><suffix>\n"
              << "\t-no_check <0/1>: disable correctness check for each blocks\n"
              << endl);
        exit(0);
    }
    InputFlags flags = InputFlags(input);
    try{
        if(flags.block_t == "fixed_size")
            dump_generic<ms_blocks>(input.getCmdOption("-ms_path"), flags);
        else if(flags.block_t == "constant_ones")
            dump_generic<const_ones_ms_blocks>(input.getCmdOption("-ms_path"), flags);
        else if(flags.block_t == "zeros_and_ones")
            dump_generic<zo_patt_ms_blocks>(input.getCmdOption("-ms_path"), flags);
        else
            throw string("Unknown block type" + flags.block_t);
    } catch (string s) {
        cerr << "ERROR:" << endl;
        cerr << s << endl;
        return 1;
    }

    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);
}


