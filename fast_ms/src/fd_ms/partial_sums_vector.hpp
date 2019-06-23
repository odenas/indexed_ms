#ifndef PARTIAL_SUMS_VECTOR_HPP
#define PARTIAL_SUMS_VECTOR_HPP

#include <bitset>
#include <string>
#include <iostream>
#include <fstream>

#include <sdsl/vectors.hpp>
#include "smsb/range_ms_sum.h"


namespace fdms {
    template<typename vec_type, typename it_type, typename size_type>
    class partial_sums_vector1 {
    public:
        const size_t m_block_size;
        const string m_path;

        partial_sums_vector1(const string path, const size_t bsize) :
            m_block_size{bsize}, m_path{path}
        { }

        static void _show_vec(vec_type& ms, it_type* it){
            cout << endl;
            for(int i=0; i<ms.getSize(); i++)
                cout << static_cast<int>(it->isSet(i)) << " ";
            cout << endl;
        }

        static size_type trivial_range_sum(vec_type& ms, it_type* it,
                size_type int_from, size_type int_to){

            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, sum_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;

            if(int_from > 0){
                bit_from = it->select(int_from - 1);
                //cout << "+ " << int_from << " -> " << bit_from << endl;
                prev_ms = bit_from - 2 * (int_from - 1);
                i = bit_from + 1;
            }
            while (cnt1 < (int_to - int_from)) {
                if (it->isSet(i)) {
                    //(cerr << "MS[" << cnt1 - 1 << "] = " << prev_ms << ", SUM = " << sum_ms << endl);
                    cur_ms = prev_ms + cnt0 - 1;
                    sum_ms += cur_ms;
                    prev_ms = cur_ms;
                    cnt0 = 0;
                    cnt1 += 1;
                } else {
                    cnt0 += 1;
                }
                i += 1;
            }
            return sum_ms;
        }

        size_type range_sum(vec_type& ms, sdsl::int_vector<64>& ridx, it_type* it,
                const size_type from, const size_type to) {

            assert(from < to);
            size_type to_sum = range_sum_prefix(ms, ridx, it, to - 1);
            size_type from_sum = (from == 0 ? 0 : range_sum_prefix(ms, ridx, it, from - 1));
            assert(from_sum <= to_sum);
            return to_sum - from_sum;
        }


        size_type range_sum_prefix(vec_type& ms, sdsl::int_vector<64>& ridx, it_type* it, const size_type to_ms_idx) {
            //cerr << "[indexed (" << flags.block_size << ")] " << to_ms_idx << endl;

            // index of last term of sum
            size_type int_ms_idx = to_ms_idx;
            size_type bit_ms_idx = it->select(int_ms_idx);
            size_type block_idx = bit_ms_idx / m_block_size;

            size_type sum_ms = 0; // to be subtracted from ridx[block_idx]
            {
                size_type prev_ms = bit_ms_idx - (2 * int_ms_idx); // needed for 1st term beyond the sum
                size_type nzeros = 0;

                // loop from bit_ms_idx + 1 to the end of the block
                for (size_type i = bit_ms_idx + 1; i < (block_idx + 1) * m_block_size; i++) {
                    if (it->isSet(i) == 1) {
                        size_type cur_ms = (prev_ms + nzeros - 1);
                        sum_ms += cur_ms; // since MS_i - MS_{i-1} + 1 = nzeros
                        prev_ms = cur_ms;
                    } else {
                        nzeros += 1;
                    }
                }
            }
            size_type answer = ridx[block_idx] - sum_ms;
            return answer;
        }

    };

    template<typename size_type, typename ms_type, typename ms_sel_1_type>
    class partial_sums_vector {
        typedef sdsl::int_vector_buffer<1> buff_vec_t;

    public:
        const size_t m_block_size;
        const string m_path;

        partial_sums_vector(const string path, const size_t bsize) :
            m_block_size{bsize}, m_path{path}
        { }

        static void _show_vec(const ms_type& ms){
            cout << endl;
            for(int i=0; i<ms.size(); i++)
                cout << ms[i] << " ";
            cout << endl;
        }

        static size_type trivial_range_sum(const ms_type& ms, ms_sel_1_type& ms_sel,
                size_type int_from, const size_type int_to) {

            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, sum_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;

            if(int_from > 0){
                bit_from = ms_sel(int_from);
                //cout << "* " << int_from << " -> " << bit_from << endl;
                prev_ms = bit_from - 2 * (int_from - 1);
                i = bit_from + 1;
            }
            while (cnt1 < (int_to - int_from)) {
                if (ms[i] == 1) {
                    //(cerr << "MS[" << cnt1 - 1 << "] = " << prev_ms << ", SUM = " << sum_ms << endl);
                    cur_ms = prev_ms + cnt0 - 1;
                    sum_ms += cur_ms;
                    prev_ms = cur_ms;
                    cnt0 = 0;
                    cnt1 += 1;
                } else {
                    cnt0 += 1;
                }
                i += 1;
            }
            return sum_ms;
        }

        static size_type djamal_range_sum(const ms_type& ms, ms_sel_1_type& ms_sel,
                size_type int_from, const size_type int_to) {

            if (int_from >= int_to)
                return 0;

            size_type bit_from = 0, bit_to = ms_sel(int_to);
            size_type prev_ms = 1;

            if(int_from > 0){
                bit_from = ms_sel(int_from);
                //cout << "* " << int_from << " -> " << bit_from << endl;
                prev_ms = bit_from - 2 * (int_from - 1);
            }
            (cerr << "prev_ms = " << prev_ms << ", "
             << "bit_from = " << bit_from << " (int_from = " << int_from << ")," 
             << "bit_to = " << bit_to << " (int_to = " << int_to << ")"
             << endl); 
            const int ss = sizeof(uint64_t);
            for(int i = 0; i <= ms.size() / ss; i+=ss){
                size_t v = ms.data()[i];
                cerr << i << ": " << std::bitset<ss>(v) << endl;
            }
            return (size_type) range_ms_sum_fast64(prev_ms, bit_from, bit_from + 1, ms.data());
        }

        size_type range_sum_prefix(ms_type& ms, sdsl::int_vector<64>& ridx, ms_sel_1_type& ms_sel, const size_type to_ms_idx) {
            //cerr << "[indexed (" << flags.block_size << ")] " << to_ms_idx << endl;

            // index of last term of sum
            size_type int_ms_idx = to_ms_idx;
            size_type bit_ms_idx = ms_sel(int_ms_idx + 1);
            size_type block_idx = bit_ms_idx / m_block_size;

            size_type sum_ms = 0; // to be subtracted from ridx[block_idx]
            {
                size_type prev_ms = bit_ms_idx - (2 * int_ms_idx); // needed for 1st term beyond the sum
                size_type nzeros = 0;

                // loop from bit_ms_idx + 1 to the end of the block
                for (size_type i = bit_ms_idx + 1; i < (block_idx + 1) * m_block_size; i++) {
                    if (ms[i] == 1) {
                        size_type cur_ms = (prev_ms + nzeros - 1);
                        sum_ms += cur_ms; // since MS_i - MS_{i-1} + 1 = nzeros
                        prev_ms = cur_ms;
                    } else {
                        nzeros += 1;
                    }
                }
            }
            size_type answer = ridx[block_idx] - sum_ms;
            return answer;
        }

        size_type range_sum(ms_type& ms, sdsl::int_vector<64>& ridx, ms_sel_1_type& ms_sel,
                const size_type from, const size_type to) {

            assert(from < to);
            size_type to_sum = range_sum_prefix(ms, ridx, ms_sel, to - 1);
            size_type from_sum = (from == 0 ? 0 : range_sum_prefix(ms, ridx, ms_sel, from - 1));
            assert(from_sum <= to_sum);
            return to_sum - from_sum;
        }

        static void dump(const string ms_path, const size_type block_size) {
            buff_vec_t ms(ms_path, std::ios::in);

            sdsl::int_vector_buffer<64> out_vec(
                    InputSpec::rdix_fname(ms_path, block_size),
                    std::ios::out);

            size_type one_cnt = 0, out_idx = 0, ms_value = 0, cum_ms = 0;
            for (size_type ms_idx = 0; ms_idx < ms.size(); ms_idx++) {
                if (ms[ms_idx] == 1) {
                    ms_value = ms_idx - 2 * one_cnt;
                    cum_ms += ms_value;
                    one_cnt += 1;
                }
                if (ms_idx and (ms_idx + 1) % block_size == 0) {
                    out_vec[out_idx++] = cum_ms;
                }
            }
        }

    };
}
#endif /* PARTIAL_SUMS_VECTOR_HPP */
