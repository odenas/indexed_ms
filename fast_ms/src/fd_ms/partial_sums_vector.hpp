#ifndef PARTIAL_SUMS_VECTOR_HPP
#define PARTIAL_SUMS_VECTOR_HPP

#include <bitset>
#include <string>
#include <iostream>
#include <fstream>

#include <sdsl/vectors.hpp>
extern "C" {
    #include "smsb/range_ms_sum.h"
    #include "smsb/naive_ms_range_sum.h"
}

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

namespace fdms {
    /* rle - based class */
    template<typename vec_type, typename it_type, typename size_type>
    class partial_sums_vector1 {
    public:

        const vec_type& m_ms;
        it_type *m_it;

        partial_sums_vector1(const vec_type& v, it_type* it) : m_ms{v}, m_it{it} {}

        void show_vec(){
            cout << endl;
            for(int i=0; i<m_ms.getSize(); i++){
                cout << (i % 10 == 0 ? "*" : " ");
            }
            cout << endl;
            for(int i=0; i<m_ms.getSize(); i++){
                cout << i % 10 << "";
            }
            cout << endl;
            for(int i=0; i<m_ms.getSize(); i++){
                cout << m_it->isSet(i) << "";
            }
            cout << endl;
        }
        size_type trivial_range_sum(const size_type int_from, const size_type int_to) {
            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, sum_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;

            if(int_from > 0){
                bit_from = m_it->select(int_from - 1);
                //cout << "+ " << int_from << " -> " << bit_from << endl;
                prev_ms = bit_from - 2 * (int_from - 1);
                i = bit_from + 1;
            }
            while (cnt1 < (int_to - int_from)) {
                if (m_it->isSet(i)) {
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
//            if(sum_ms != rle_range_sum(int_from, int_to)){
//                cerr << "rle_range_sum(" << int_from << ", " << int_to << ") != " << sum_ms << endl;
//                exit(1);
//            }
            return sum_ms;
        }

        size_type rle_range_sum(const size_type int_from, const size_type int_to){
            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, sum_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;

            if(int_from > 0){
                bit_from = m_it->select(int_from - 1);
                //cout << "+ " << int_from << " -> " << bit_from << endl;
                prev_ms = bit_from - 2 * (int_from - 1);
                i = bit_from;
            }
            // (i, l). i = index of next 1, l = nr. of remaining 1s in that run
            // std::pair<size_type, size_type> it->selectNextRun(m_ms.getSize());

            while (cnt1 < (int_to - int_from)) {
                size_type j = m_it->selectNext();
                if(j <= i) {
                    show_vec();
                    cerr << "[" << m_ms.getSize() << "] i = " << i << ", but j = " << j << endl;
                    exit(1);
                }
                cnt0 = j - i - 1;
                cur_ms = prev_ms + cnt0 - 1;
                sum_ms += cur_ms;
                prev_ms = cur_ms;
                cnt1 += 1;
                i = j;
            }
            return sum_ms;
        }

        size_type range_sum(sdsl::int_vector<64>& ridx,
                const size_type from, const size_type to, const size_type bsize) {

            assert(from < to);
            size_type to_sum = range_sum_prefix(ridx, to - 1, bsize);
            size_type from_sum = (from == 0 ? 0 : range_sum_prefix(ridx, from - 1, bsize));
            assert(from_sum <= to_sum);
            return to_sum - from_sum;
        }

        size_type range_sum_prefix(sdsl::int_vector<64>& ridx, const size_type to_ms_idx, const size_type bsize) {
            //cerr << "[indexed (" << flags.block_size << ")] " << to_ms_idx << endl;

            // index of last term of sum
            size_type int_ms_idx = to_ms_idx;
            size_type bit_ms_idx = m_it->select(int_ms_idx);
            size_type block_idx = bit_ms_idx / bsize;
            size_type sum_ms = 0; // to be subtracted from ridx[block_idx]
            {
                size_type prev_ms = bit_ms_idx - (2 * int_ms_idx); // needed for 1st term beyond the sum
                size_type nzeros = 0;

                // loop from bit_ms_idx + 1 to the end of the block
                for (size_type i = bit_ms_idx + 1; i < (block_idx + 1) * bsize; i++) {
                    if (m_it->isSet(i) == 1) {
                        size_type cur_ms = (prev_ms + nzeros - 1);
                        sum_ms += cur_ms; // since MS_i - MS_{i-1} + 1 = nzeros
                        prev_ms = cur_ms;
                        nzeros = 0;
                    } else {
                        nzeros += 1;
                    }
                }
            }
            size_type answer = ridx[block_idx] - sum_ms;
            return answer;
        }
    };

    /* sdsl based class */
    template<typename vec_type, typename ms_sel_1_type, typename size_type>
    class partial_sums_vector {
        typedef sdsl::int_vector_buffer<1> buff_vec_t;

    public:
        const vec_type& m_ms;
        ms_sel_1_type& m_ms_sel;

        partial_sums_vector(const vec_type& ms, ms_sel_1_type& ms_sel) :
            m_ms{ms}, m_ms_sel{ms_sel} { }
        void _show_vec(){
            cout << endl;
            for(int i=0; i<m_ms.size(); i++)
                cout << m_ms[i] << " ";
            cout << endl;
        }

        size_type trivial_range_sum(const size_type int_from, const size_type int_to) {
            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, sum_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;

            if(int_from > 0){
                bit_from = m_ms_sel(int_from);
                //cout << "* " << int_from << " -> " << bit_from << endl;
                prev_ms = bit_from - 2 * (int_from - 1);
                i = bit_from + 1;
            }
            while (cnt1 < (int_to - int_from)) {
                if (m_ms[i] == 1) {
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

        /**
         * Sum the values MS[i] for int_from <= i < int_to
         */
        size_type djamal_range_sum(const size_type int_from, const size_type int_to) {
            if (int_from >= int_to)
                return 0;

            size_type bit_from = m_ms_sel(int_from + 1);
            size_type bit_to = m_ms_sel(int_to);
            size_type prev_ms = 1;

            //cout << "* " << int_from << " -> " << bit_from << endl;
            prev_ms = bit_from - 2 * int_from;
            //(cerr << "prev_ms = " << prev_ms << ", "
            // << "bit_from = " << bit_from << " (int_from = " << int_from << "),"
            // << "bit_to = " << bit_to << " (int_to = " << int_to << ")"
            // << endl);
            const int ss = sizeof(size_t);
            return (size_type) range_ms_sum_fast64(prev_ms, bit_from, bit_to, m_ms.data());
            //return (size_type) naive_range_ms64(int_from, int_to - 1, 2048, ms.data());
        }

        size_type range_sum_prefix(sdsl::int_vector<64>& ridx, const size_type to_ms_idx, const size_type bsize) {
            //cerr << "[indexed (" << flags.block_size << ")] " << to_ms_idx << endl;

            // index of last term of sum
            size_type int_ms_idx = to_ms_idx;
            size_type bit_ms_idx = m_ms_sel(int_ms_idx + 1);
            size_type block_idx = bit_ms_idx / bsize;

            size_type sum_ms = 0; // to be subtracted from ridx[block_idx]
            {
                size_type prev_ms = bit_ms_idx - (2 * int_ms_idx); // needed for 1st term beyond the sum
                size_type nzeros = 0;

                // loop from bit_ms_idx + 1 to the end of the block
                for (size_type i = bit_ms_idx + 1; i < (block_idx + 1) * bsize; i++) {
                    if (m_ms[i] == 1) {
                        size_type cur_ms = (prev_ms + nzeros - 1);
                        sum_ms += cur_ms; // since MS_i - MS_{i-1} + 1 = nzeros
                        prev_ms = cur_ms;
                        nzeros = 0;
                    } else {
                        nzeros += 1;
                    }
                }
            }
            size_type answer = ridx[block_idx] - sum_ms;
            return answer;
        }

        size_type range_sum(sdsl::int_vector<64>& ridx,  const size_type from, const size_type to, const size_type bsize) {
            assert(from < to);
            size_type to_sum = range_sum_prefix(ridx, to - 1, bsize);
            size_type from_sum = (from == 0 ? 0 : range_sum_prefix(ridx, from - 1, bsize));
            assert(from_sum <= to_sum);
            return to_sum - from_sum;
        }

        static void dump(const string ms_path, const size_type block_size) {
            buff_vec_t ms(ms_path, std::ios::in);
            if(ms.size() % block_size != 0)
                throw string{"block_size (" + to_string(block_size) + ")" +
                             " should divide the size of the given ms " +
                             "(" + to_string(ms.size()) +")"};

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
                    cerr << cum_ms << " ";
                }
                if(ms_idx == 125)
                    cerr << ".";
            }
            cerr << endl;
        }

    };
}
#endif /* PARTIAL_SUMS_VECTOR_HPP */
