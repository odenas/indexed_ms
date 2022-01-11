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
    bool check;

    InputFlags() { }

    InputFlags(const InputFlags& f) :
        start{f.start}, len{f.len}
    { }

    InputFlags(const size_type start, const size_type len, const bool int_format) :
        start{start}, len{len}
    { }

    InputFlags(OptParser input) :
        start{static_cast<size_type> (std::stoll(input.getCmdOption("-start")))},
        len{static_cast<size_type> (std::stoll(input.getCmdOption("-len")))},
        check{static_cast<bool> (std::stoll(input.getCmdOption("-check")))} {}
};

void sel_result_is_valid(sdsl::bit_vector& ms, size_type i, size_type sel_i){
}

size_type _select(sdsl::bit_vector& ms, sdsl::bit_vector::select_1_type& ms_sel, const size_type i){
    size_type sel_i = ms_sel(i + 1);  // sdsl-select is 1-based
    if(ms[sel_i] == 0)
        throw string{"ms["} + to_string(sel_i) + string{"] = 0"};
    return sel_i;
}

size_type _int_ms(size_type sel_i, size_type i){
    if(sel_i < 2 * i)
        throw (string{"i = "} + to_string(i) +
               string{" sel(i) = "} + to_string(sel_i) +
               string{", but 2i = "} + to_string(2 * i));
    return sel_i - 2 * i;
}

int check(const string ms_path){
    cout << "checking: " << ms_path << endl;
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);
    size_type query_size = ms.size() / 2;
    sdsl::bit_vector::select_1_type ms_sel(&ms);

    size_type prev_ms = 1;
    size_type prev_bit_i = 0;
    try{
        for(size_type i = 0; i < query_size; i++){
            if(i > 0){
                // compute using the formula MS[i] = sel(i) - 2i
                prev_bit_i = _select(ms, ms_sel, i - 1);
                prev_ms = _int_ms(prev_bit_i, i - 1);
            }

            // compute using the formula MS[i] = sel(i) - 2i
            size_type curr_bit_i = _select(ms, ms_sel, i);
            size_type curr_ms = _int_ms(curr_bit_i, i);

            if (i > 0 and curr_bit_i <= prev_bit_i){
                throw (string{"i = "} + to_string(i) +
                       string{" sel(i) = "} + to_string(curr_ms) +
                       string{", but sel(i - 1) = "} + to_string(prev_ms));
            }

            size_type n_zeros = curr_bit_i - prev_bit_i - (i > 0);

            // check that MS[i] - MS[i-1] + 1 = # zeros between the two corresponding 1s in ms
            if(curr_ms + 1 - prev_ms != n_zeros){
                throw (string{"i = "} + to_string(i) +
                       string{" curr_ms = "} + to_string(curr_ms) +
                       string{", prev_ms = "} + to_string(prev_ms) +
                       string{" but #0s between the two = "} + to_string(curr_bit_i - prev_bit_i - (prev_bit_i > 0)));
            }
        }
    } catch (string s){
        throw s;
    }
    return 0;
}

int comp(const string ms_path, const InputFlags& flags) {
    buff_vec_t ms(ms_path, std::ios::in);

    size_type end = flags.start + flags.len;
    if(flags.len == 0)
        end = ms.size();

    size_type sum = 0;

    for(size_type j = flags.start; j < end; j++){
        sum += ms[j];
        if(j >= ms.size()){
            cout << "reached the end at " << j << endl;
            break;
        }
    }
    cout << "filename:    " << ms_path << endl;
    cout << "total size:  " << ms.size() << endl;
    cout << "range:       [" << flags.start << ", " << end << ")" << endl;
    cout << "range size:  " << end - flags.start << endl;
    cout << "0s in range: " << end - flags.start - sum << endl;
    cout << "1s in range: " << sum << endl;
    if(flags.check){
        try{
            return check(ms_path);
        } catch (string s){
            throw s;
        }
    }

    return 0;
}

int main(int argc, char **argv) {
    if(argc == 1){
        (cerr << "Print a section of the ms bit-vector\n"
              << "Args:\n"
              << help__ms_path
              << "\t-start <non-negative int>: start of a 0-based half-open interval [from_idx, to_idx)\n"
              << "\t-len <non-negative int>: length of a 0-based half-open interval [from_idx, to_idx)\n"
              << "\t-check <0/1>: check consistency\n"
              << endl);
        exit(0);
    }

    OptParser input(argc, argv);
    string ms_path = input.getCmdOption("-ms_path");

    InputFlags flags;
    try{
        flags = InputFlags(input);
        return comp(ms_path, flags);
    } catch (string s) {
        cerr << s << endl;
        return 1;
    }
}

