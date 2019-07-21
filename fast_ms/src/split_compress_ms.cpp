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

#include "../malloc_count/malloc_count.h"


using namespace fdms;
using namespace std;

typedef typename ms_blocks::size_type size_type;
typedef typename ms_blocks::pair_t pair_t;


class Counter{
public:
    map<string, size_type> reg;

    string ms_block_raw = "ms_block";
    string ms_block_compr = "ms_block_compressed";
    string ms_block_check = "ms_block_check";
    string select = "select1";
    string ms_total = "ms_total";

    Counter(){
        reg[ms_block_check] = 0;
    }

    size_type operator[](string key) {
        return reg[key];
    }

    size_type abs_point() const {
        malloc_count_reset_peak();
        return (size_type) malloc_count_peak();
    }

    void block_reset(){
        vector<string> keys = {ms_block_raw, ms_block_check, ms_block_compr, select};
        for(auto key: keys)
            reg[key] = 0;
    }

    size_type add(const string key, const size_type from){
        size_type to = abs_point();
        if (from > to)
            throw string{"from peak (" + to_string(from) + ") < to peak (" + to_string(to) + ")"};
        reg[key] = to - from;
        return reg[key];
    }
};


class InputFlags {
public:
    size_type len, cnt1;
    string block_t;
    bool check;

    InputFlags() { }

    InputFlags(const InputFlags& f) : cnt1{f.cnt1}, len{f.len}, check{f.check} { }

    InputFlags(OptParser input) :
            len{static_cast<size_type> (std::stoll(input.getCmdOption("-len")))},
            cnt1{static_cast<size_type> (std::stoll(input.getCmdOption("-cnt1")))},
            check{input.getCmdOption("-no_check") == "0"}
    {
        if(len == 0 || cnt1 == 0){
            throw string{"needs -len and -cnt1"};
        }
    }

    size_type block_l() const {
        if(block_t == "constant_ones")
            return cnt1;
        else if(block_t == "fixed_size")
            return len;
        throw string{"Bad block type: '" + block_t + "'."};
    }
};


template<typename enc_type>
size_type fill_encoder(const sdsl::bit_vector& ms, enc_type& encoder){
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

template<typename block_type_cls>
void subvector_footprint(const sdsl::bit_vector& ms,
        const size_type slice_idx, block_type_cls& blocks,
        const string& compr,
        const bool check,
        Counter& mem_usage){

    size_type start = blocks.slices[slice_idx].first;
    size_type end = blocks.slices[slice_idx].second;
    size_type abs_point = 0;

    if(start >= ms.size())
        throw string{"Out of bounds on ms vector (len " + std::to_string(ms.size()) + "). " +
                     "start=" + std::to_string(start)};

    //cerr << " *** [" << start << ", " << end << ")" << endl;

    abs_point = mem_usage.abs_point();
    size_type w = 64, len = end - start;
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
    mem_usage.add(mem_usage.ms_block_raw, abs_point);

    if(check){
        abs_point = mem_usage.abs_point();
        cerr << "checking [" << start << ", " << end << ") ..." << endl;
        string out_path = compr + ".ms";
        sdsl::store_to_file(out_ms, out_path);
        try{
            blocks.check(out_path, slice_idx);
        } catch (string s) {
            throw s;
        }
        mem_usage.add(mem_usage.ms_block_check, abs_point);
    }

    string out_path = compr + ".ms";
    if (compr == "baseline"){
        // no compression
        abs_point = mem_usage.abs_point();
        mem_usage.add(mem_usage.ms_block_compr, abs_point);

        // no select1
        abs_point = mem_usage.abs_point();
        mem_usage.add(mem_usage.select, abs_point);

        //save
        //sdsl::store_to_file(out_ms, out_path);
    }
    else if(compr == "no_compr"){
        // no compression
        abs_point = mem_usage.abs_point();
        mem_usage.add(mem_usage.ms_block_compr, abs_point);

        // select1
        abs_point = mem_usage.abs_point();
        sdsl::bit_vector::select_1_type a(&out_ms);
        mem_usage.add(mem_usage.select, abs_point);

        //save
        //sdsl::store_to_file(out_ms, out_path);
    }
    else if(compr == "rrr"){
        //compression
        abs_point = mem_usage.abs_point();
        sdsl::rrr_vector<> compr_out_ms(out_ms);
        compr_out_ms.size();
        mem_usage.add(mem_usage.ms_block_compr, abs_point);

        // select
        abs_point = mem_usage.abs_point();
        sdsl::rrr_vector<>::select_1_type a(&compr_out_ms);
        a(1);
        mem_usage.add(mem_usage.select, abs_point);

        //save
        //sdsl::store_to_file(compr_out_ms, out_path);
    }
    else if (compr == "rle"){
        //compression
        abs_point = mem_usage.abs_point();
        CSA::RLEEncoder encoder(32);
        size_type n_runs = fill_encoder<CSA::RLEEncoder>(out_ms, encoder);
        mem_usage.add("rle_encoder", abs_point);
        abs_point = mem_usage.abs_point();
        CSA::RLEVector compr_out_ms(encoder, ms.size());
        mem_usage.add(mem_usage.ms_block_compr, abs_point);
        mem_usage.reg["n_runs"] = n_runs;

        // no select
        abs_point = mem_usage.abs_point();
        mem_usage.add(mem_usage.select, abs_point);

        //save
        //std::ofstream out{out_path, std::ios::binary}; compr_out_ms.writeTo(out);
    }
    else
        throw string{"Bad compression string: " + compr};
}

template<typename block_type_cls>
int dump_generic(const string ms_path, const InputFlags& flags){
    cerr << flags.block_t << endl;

    vector<string> compr_types = {"baseline", "no_compr", "rle", "rrr"};
    for(auto compr: compr_types){
        cerr << " * " << compr << "..." << endl;
        Counter mem_usage{};

        // all of ms footprint
        size_type abs_point = mem_usage.abs_point();
        sdsl::bit_vector ms;
        sdsl::load_from_file(ms, ms_path);
        mem_usage.add(mem_usage.ms_total, abs_point);
        cerr << " ** loaded input of size " << ms.size() << endl;

        block_type_cls blocks(ms, flags.block_l());
        cerr << " ** computed " << blocks.slices.size() << " blocks (" << flags.block_t << ")" << endl;
        for(size_type slice_idx = 0; slice_idx < blocks.slices.size(); slice_idx++){
            //cerr << " *** " << slice_idx << endl;

            subvector_footprint<block_type_cls>(ms, slice_idx, blocks, compr, flags.check, mem_usage);
            for (auto item : mem_usage.reg){
                (cout << flags.block_t <<
                 "," << blocks.slices[slice_idx].first <<
                 "," << blocks.slices[slice_idx].second <<
                 "," << compr <<
                 "," << item.first << "," << item.second << endl);
            }
            mem_usage.block_reset();
        }
        cerr << " * DONE" << endl;
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
              << endl);
        exit(0);
    }
    try{
        InputFlags flags = InputFlags(input);
        cout << "block_type,start,end,compr,mempoint,value_bytes" << endl;

        flags.block_t = "fixed_size";
        dump_generic<ms_blocks>(input.getCmdOption("-ms_path"), flags);

        flags.block_t = "constant_ones";
        dump_generic<const_ones_ms_blocks>(input.getCmdOption("-ms_path"), flags);
    } catch (string s) {
        cerr << "ERROR: " << endl << s << endl;
        return 1;
    }
    return 0;
}
