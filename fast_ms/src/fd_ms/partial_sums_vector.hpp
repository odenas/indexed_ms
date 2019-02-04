/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   partial_sums_vector.hpp
 * Author: brt
 *
 * Created on February 3, 2019, 9:57 PM
 */

#ifndef PARTIAL_SUMS_VECTOR_HPP
#define PARTIAL_SUMS_VECTOR_HPP

#include <string>
#include <iostream>
#include <fstream>

#include <sdsl/vectors.hpp>


namespace fdms {

    template<typename size_type>
    class partial_sums_vector {
        typedef sdsl::int_vector_buffer<1> buff_vec_t;

    public:
        const size_t m_block_size;
        const string m_path;

        partial_sums_vector(const string path, const size_t bsize) : m_block_size{bsize}, m_path{path}
        {
        }

        size_type range_sum(sdsl::bit_vector& ms, sdsl::int_vector<64>& ridx,
                sdsl::bit_vector::select_1_type ms_sel, sdsl::bit_vector::rank_1_type ms_rank,
                const size_type from, const size_type to) {

            assert(from < to);
            size_type to_sum = range_sum_prefix(ms, ridx, ms_sel, ms_rank, to - 1);
            size_type from_sum = range_sum_prefix(ms, ridx, ms_sel, ms_rank, from);
            assert(from_sum <= to_sum);
            return to_sum - from_sum;
        }

        size_type range_sum_prefix(sdsl::bit_vector& ms, sdsl::int_vector<64>& ridx,
                sdsl::bit_vector::select_1_type ms_sel, sdsl::bit_vector::rank_1_type ms_rank,
                const size_type to_ms_idx) {
            //cerr << "[indexed (" << flags.block_size << ")] " << to_ms_idx << endl;

            // index of last term of sum
            size_type int_ms_idx = to_ms_idx;
            size_type bit_ms_idx = ms_sel(int_ms_idx + 1);
            size_type block_idx = bit_ms_idx / m_block_size;

            size_type sum_ms = 0; // to be subtracted from ridx[block_idx]
            {
                size_type prev_ms = bit_ms_idx - (2 * int_ms_idx); // needed for 1st term beyond the sum
                size_type end_block_nones = ms_rank((block_idx + 1) * m_block_size);
                size_type nzeros = 0, nones = int_ms_idx;

                // loop from bit_ms_idx + 1 to the end of the block
                for (size_type i = bit_ms_idx + 1; i < (block_idx + 1) * m_block_size; i++) {
                    if (nones >= end_block_nones)
                        break;

                    if (ms[i] == 1) {
                        size_type cur_ms = (prev_ms + nzeros - 1);
                        sum_ms += cur_ms; // since MS_i - MS_{i-1} + 1 = nzeros
                        nones += 1;
                        prev_ms = cur_ms;
                    } else {
                        nzeros += 1;
                    }
                }
            }
            size_type answer = ridx[block_idx] - sum_ms;
            return answer;
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

        static size_type trivial_range_sum(const sdsl::bit_vector& ms, sdsl::bit_vector::select_1_type ms_sel,
                size_type int_from, const size_type int_to) {

            size_type bit_from = ms_sel(int_from + 1);
            size_type prev_ms = 1, cur_ms = 0, sum_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;
            while (cnt1 < (int_to - int_from)) {
                if (ms[i] == 1) {
                    //(cout << "MS[" << cnt1 - 1 << "] = " << prev_ms << ", SUM = " << sum_ms << endl);
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
    };
}
#endif /* PARTIAL_SUMS_VECTOR_HPP */

