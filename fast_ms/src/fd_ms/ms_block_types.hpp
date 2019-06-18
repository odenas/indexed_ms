//
//  ms_block_types.hpp
//  fast_ms
//
//  Created by denas on 10/27/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef ms_block_types_h
#define ms_block_types_h

#include <cassert>
#include <vector>
#include <string>
#include "sdsl/bit_vectors.hpp"


using namespace std;


namespace fdms {
    class ms_blocks {
    public:
        typedef long unsigned int size_type;
        typedef pair<size_type, size_type> pair_t;

        vector<pair_t> slices;

        ms_blocks() = default;

        size_type size(size_type i) const {
            return slices[i].second - slices[i].first;
        }

        ms_blocks(const sdsl::bit_vector& ms, const size_type slice_size){
            size_type input_size = ms.size();
            size_type chunks = input_size / slice_size;
            size_type extra = input_size % slice_size;
            slices = vector<pair_t>(chunks + (extra != 0));
            for (size_type i = 0, from = 0; i < chunks; i++, from+=slice_size) {
                slices[i] = make_pair(from, from + slice_size);
            }
            if(slices.size() > chunks){
                assert(extra != 0);
                slices[chunks] = make_pair(input_size - extra, input_size);
            }
        }

        string make_out_fname(const string prefix, const string suffix, const size_type idx) const {
            size_type start = slices[idx].first;
            size_type end = slices[idx].second;
            return (prefix + "_" + to_string(start) + "_" + to_string(end) + suffix);
        }

        void check(const string ms_path, const size_type slice_idx) const {
            sdsl::bit_vector block_ms;
            sdsl::load_from_file(block_ms, ms_path);
            try{check(block_ms, slice_idx);}
            catch (string s){throw s;}
        }

        void check(const sdsl::bit_vector& ms, const size_type slice_idx) const {
            //cout << "[" << slices[slice_idx].first << ", " << slices[slice_idx].second << ")" << endl;
            size_type expected_len = slices[slice_idx].second - slices[slice_idx].first;
            if(ms.size() != expected_len)
                throw string{"Expecting length " + to_string(expected_len) +
                             ". Was " + to_string(ms.size()) + "."};
        }
    };

    class const_ones_ms_blocks : public ms_blocks {
    public:
        vector<size_type> n_ones;

        const_ones_ms_blocks() {
            slices = vector<pair_t>(0);
            n_ones = vector<size_type>(0);
        }

        const_ones_ms_blocks(const sdsl::bit_vector& ms, const size_type len) {
            slices = vector<pair_t>(0);
            n_ones = vector<size_type>(0);

            size_type start = 0, end = 0;
            size_type cnt1 = 0;
            do{
                if(ms[end++] == 1)
                    cnt1 += 1;

                if(cnt1 == len){
                    slices.push_back(make_pair(start, end));
                    n_ones.push_back(cnt1);
                    start = end;
                    cnt1 = 0;
                }
            } while(end < ms.size());
        }

        void check(const sdsl::bit_vector& ms, const size_type slice_idx) const {
            size_type expected_len = slices[slice_idx].second - slices[slice_idx].first;
            size_type expected_n_ones = n_ones[slice_idx];

            if(ms.size() != expected_len)
                throw string{"Expecting length " + to_string(expected_len) +
                             ". Was " + to_string(ms.size()) + "."};
            size_type c = 0;
            for(size_type i = 0; i < ms.size(); i++)
                c += ms[i];
            if(c < expected_n_ones)
                throw string{"Expecting " + to_string(expected_len) +
                             " ones. Found " + to_string(c) + "."};
        }
        using ms_blocks::check;
    };

    class zo_patt_ms_blocks : public const_ones_ms_blocks {
        typedef sdsl::bit_vector::select_0_type ms_sel_0_type;
        typedef sdsl::bit_vector::select_1_type ms_sel_1_type;
    public:
        zo_patt_ms_blocks(const sdsl::bit_vector& ms, const size_type min_len){
            const_ones_ms_blocks();

            ms_sel_0_type ms_sel0(&ms);
            ms_sel_1_type ms_sel1(&ms);

            size_type cnt0 = 1, cnt1 = 1;  // initialize cummulative cnt of 0s and 1s
            size_type idx0 = ms_sel0(1), idx1 = ms_sel1(1);  // initialize indexes of 0s and 1s

            assert(idx1 != idx0);
            if(idx1 < idx0){  // this vector started with 1s
                assert(idx1 == 0);
                cnt1 += (idx0 - idx1);
                idx1 = ms_sel1(cnt1);
            }
            assert(idx1 > idx0);
            //cout << cnt1 << "; " << idx1 << endl;

            do{
                size_type c = ms_sel0(cnt0 + idx1 - idx0);
                if(c - idx1 > ms.size()) // overflowed
                    break;
                //cout << "c: " << c << endl;
                if(c - idx0 >= min_len){
                    slices.push_back(make_pair(idx0, c));
                    n_ones.push_back(c - idx1);
                }
                cnt0 += idx1 - idx0;
                idx0 = c;
                //cout << cnt0 << "; " << idx0 << endl;
                cnt1 += c - idx1;
                idx1 = ms_sel1(cnt1);
            } while(cnt0 + cnt1 < ms.size());
        }

        void check(const sdsl::bit_vector& ms, const size_type slice_idx) const {
            size_type expected_len = slices[slice_idx].second - slices[slice_idx].first;
            size_type expected_n_ones = n_ones[slice_idx];

            if(ms.size() != expected_len)
                throw string{"Expecting length " + to_string(expected_len) +
                             ". Was " + to_string(ms.size()) + "."};

            ms_sel_1_type ms_sel1(&ms);
            size_type idx1 = ms_sel1(1);
            for(size_type i = 0; i < idx1; i++){  // expecting a string of 0s
                if(ms[i] != 0)
                    throw string{"Expecting a string of 0s from 0 to " +
                                 to_string(idx1) + ". v[" + to_string(i) + "] = 1."};
            }

            for(size_type i = idx1; i < ms.size(); i++){  // expecting a string of 1s
                if(ms[i] != 1)
                    throw string{"Expecting a string of 1s from " +
                                 to_string(idx1) + " to " +
                                 to_string(ms.size()) + ". v[" + to_string(i) + "] = 0."};
            }
        }
        using ms_blocks::check;
    };
};

#endif /* ms_block_types */

