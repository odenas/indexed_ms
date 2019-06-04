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
#include "fd_ms/p_ms_vector.hpp"
#include "fd_ms/stree_sct3.hpp"


using namespace fdms;
using namespace std;

typedef StreeOhleb<> cst_t;
typedef typename cst_t::size_type size_type;
typedef sdsl::bit_vector::select_0_type ms_sel_0_type;
typedef sdsl::bit_vector::select_1_type ms_sel_1_type;


class ms_blocks {
public:
    enum block_types {
        fixed_size, constant_ones, zeros_and_ones
    };

private:
    static std::map<block_types, string> c2s() {
        std::map<block_types, string> c2s = {
            {block_types::fixed_size, ".fixed_size"},
            {block_types::constant_ones, ".constant_ones"},
            {block_types::zeros_and_ones, ".zeros_and_ones"},
        };
        return c2s;
    }

public:
    ms_blocks() = default;

    static string to_str(const block_types ct){
        return c2s()[ct];
    }

    static block_types parse_block_type(const string c_str){
        for(auto item: c2s()){
            if(item.second == ("." + c_str))
                return item.first;
        }
        throw string{"bad compression string: " + c_str};
    }
};


class InputFlags {
public:
    ms_blocks::block_types blocks;
    size_type len;
    bool check;

    InputFlags() { }

    InputFlags(const InputFlags& f) :
            blocks{f.blocks}, len{f.len}, check{f.check} { }

    InputFlags(const ms_blocks::block_types blocks, const size_type len, const bool check) :
            blocks{blocks}, len{len}, check{check} { }

    InputFlags(OptParser input) :
            len{static_cast<size_type> (std::stoll(input.getCmdOption("-len")))},
            check{input.getCmdOption("-no_check") == "0"}
    {
        blocks = ms_blocks::parse_block_type(input.getCmdOption("-block_type"));
    }
};

void dump_subvector(const sdsl::bit_vector& ms, const size_type start, const size_type end, const string out_path, const bool check){
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

    if(check){
        for(size_type i = start; i < end; i++){
            if(out_ms[i - start] != ms[i])
                throw string{"Check error"};
        }
    }
    sdsl::store_to_file(out_ms, out_path);
}

string make_out_fname(const string prefix, const string suffix, const size_type start, const size_type size){
    return (prefix + "_" + to_string(start) + "_" + to_string(size) + suffix);
}

void check_fixed_size(const string path, const size_type s){
    sdsl::bit_vector block_ms;
    sdsl::load_from_file(block_ms, path);

    if(s != block_ms.size())
        throw string{"Bad size"};
}
int dump_fixed_size(const sdsl::bit_vector& ms, const string prefix, const string suffix, const InputFlags flags){
    cerr << ms.size() << "..." << flags.len << endl;
    BlockSlices<size_type> slices(ms.size(), flags.len);
    for(int slice_idx = 0; slice_idx < slices.slices.size(); slice_idx++){
        string out_path = make_out_fname(prefix, suffix, slices[slice_idx].first, slices[slice_idx].second);
        dump_subvector(ms, slices[slice_idx].first, slices[slice_idx].second, out_path, flags.check);
        if(flags.check){
            try{
                check_fixed_size(out_path, slices[slice_idx].second - slices[slice_idx].first);
            } catch (string s){
                throw s;
            }
        }
    }
    return 0;
}

void check_constant_ones(const string path, const size_type s, const size_type ones){
    sdsl::bit_vector block_ms;
    sdsl::load_from_file(block_ms, path);

    if(s != block_ms.size())
        throw string{"Bad size"};

    size_type cnt = 0;
    for(size_type i = 0; i < block_ms.size(); i++)
        cnt += block_ms[i];
    if(cnt != ones)
        throw string{"expecting a different number of ones"};
}
int dump_constant_ones(const sdsl::bit_vector& ms, const string prefix, const string suffix, const InputFlags flags){
    size_type start = 0, end = 0;
    size_type cnt1 = 0;
    do{
        if(ms[end++] == 1)
            cnt1 += 1;

        if(cnt1 == flags.len){
            string out_path = make_out_fname(prefix, suffix, start, end - start);
            dump_subvector(ms, start, end, out_path, flags.check);
            if(flags.check){
                try{
                    check_constant_ones(out_path, start - end, cnt1);
                } catch (string s){
                    throw s;
                }
            }
            start = end;
            cnt1 = 0;
        }
    }while(end < ms.size());
    return 0;
}

void check_zeros_and_ones(const string path, const size_type s){
    sdsl::bit_vector block_ms;
    sdsl::load_from_file(block_ms, path);

    if(s != block_ms.size())
        throw string{"Bad size"};

    ms_sel_1_type ms_sel1(&block_ms);
    for(size_type i = 0; i < ms_sel1(1); i++){
        if(block_ms[i] != 0)
            throw string{"expecting 0s in this part"};
    }
    for(size_type i = ms_sel1(1); i < block_ms.size(); i++){
        if(block_ms[i] != 1)
            throw string{"expecting 1s in this part"};
    }
}
int dump_zeros_and_ones(const sdsl::bit_vector& ms, const string prefix, const string suffix, const InputFlags flags){
    ms_sel_0_type ms_sel0(&ms);
    ms_sel_1_type ms_sel1(&ms);

    size_type cnt0 = 1;
    size_type idx0 = ms_sel0(cnt0);
    //cout << cnt0 << "; " << idx0 << endl;

    size_type cnt1 = 1;
    size_type idx1 = ms_sel1(1);
    if(idx1 < idx0)
        idx1 = ms_sel1(++cnt1);
    //cout << cnt1 << "; " << idx1 << endl;

    size_type c = 0;
    do{
        c = ms_sel0(cnt0 + idx1 - idx0);
        if(c - idx1 > ms.size()) // overflowed
            break;
        //cout << "c: " << c << endl;
        if(c - idx0 >= flags.len){
            string out_path = make_out_fname(prefix, suffix, idx0, c-idx0);
            dump_subvector(ms, idx0, c - idx0, out_path, flags.check);
            if(flags.check){
                try{
                    check_zeros_and_ones(out_path, c - idx0);
                } catch (string s){
                    throw s;
                }
            }
        }
        /*
        for(int i = idx0; i < c; i++){
            cout << i << " : " << ms[i] << endl;
        }
        cout << "--" << endl;
        for(int i = c; i < c + 2; i++){
            cout << i << " : " << ms[i] << endl;
        }
        cout << endl << endl;
        */

        cnt0 += idx1 - idx0;
        idx0 = ms_sel0(cnt0);
        //cout << cnt0 << "; " << idx0 << endl;
        cnt1 += c - idx1;
        idx1 = ms_sel1(cnt1);
        //cout << cnt1 << "; " << idx1 << endl;
        //if(c % 1000000 == 0)
        //    cerr << 100 * c / ms.size() << "%" << endl;
    }while(cnt0 + cnt1 < ms.size());
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
    InputFlags flags;
    try{
        flags = InputFlags(input);
    }
    catch (string s) {
        cerr << s << endl;
        return 1;
    }
    try{
        string prefix = input.getCmdOption("-out_prefix");
        string suffix = input.getCmdOption("-out_suffix");
        sdsl::bit_vector ms;
        sdsl::load_from_file(ms, input.getCmdOption("-ms_path"));

        switch(flags.blocks)
        {
            case ms_blocks::block_types::fixed_size:
                return dump_fixed_size(ms, prefix, suffix, flags);
            case ms_blocks::block_types::constant_ones:
                return dump_constant_ones(ms, prefix, suffix, flags);
            case ms_blocks::block_types::zeros_and_ones:
                return dump_zeros_and_ones(ms, prefix, suffix, flags);
            default:
                cerr << "Error." << endl;
                return 1;
        }

    } catch (string s) {
        cerr << "ERROR from dump_subvector():" << endl;
        cerr << s << endl;
    }
}

