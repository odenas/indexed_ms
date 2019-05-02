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
    ms_sel_0_type ms_sel0(&ms);
    ms_sel_1_type ms_sel1(&ms);
    cerr << "DONE" << endl;

    size_type cnt0 = 1;
    size_type idx0 = ms_sel0(cnt0);
    //cout << cnt0 << "; " << idx0 << endl;

    size_type cnt1 = 1;
    size_type idx1 = ms_sel1(1);
    if(idx1 < idx0)
        idx1 = ms_sel1(++cnt1);
    //cout << cnt1 << "; " << idx1 << endl;

    cout << "start,cnt0,cnt1,cnt" << endl;
    size_type c = 0;
    do{
        c = ms_sel0(cnt0 + idx1 - idx0);
        if(c - idx1 > ms.size()) // overflowed
            break;
        //cout << "c: " << c << endl;
        if(c - idx0 >= flags.len){
            cout << idx0 << "," << idx1 - idx0 << "," << c - idx1 << "," << c - idx0 << endl;
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
        if(c % 1000000 == 0)
            cerr << 100 * c / ms.size() << "%" << endl;
    }while(cnt0 + cnt1 < ms.size());
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