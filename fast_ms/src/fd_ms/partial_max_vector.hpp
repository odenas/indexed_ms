#ifndef PARTIAL_MAX_VECTOR_HPP
#define PARTIAL_MAX_VECTOR_HPP

#include <bitset>
#include <string>
#include <iostream>
#include <fstream>

#include <sdsl/vectors.hpp>
#include <sdsl/rmq_support.hpp>
extern "C" {
    #include "virtual_smsb/ms_range_max.h"
    #include "virtual_smsb/naive_ms_range_max.h"
}
#include "counter.hpp"
#include "range_query.hpp"

namespace fdms {

    /* sdsl based class */
    template<typename vec_type, typename ms_sel_1_type, typename size_type>
    class sdsl_partial_max_vector {
    protected:
        typedef sdsl::int_vector_buffer<1> buff_vec_t;
        typedef Counter<size_type> counter_t;
    public:
        vec_type m_ms;
        ms_sel_1_type m_ms_sel;

        sdsl_partial_max_vector(const string& ms_path) {
            sdsl::load_from_file(m_ms, ms_path);
            m_ms_sel = ms_sel_1_type(&m_ms);
        }

        sdsl_partial_max_vector(const string& ms_path, counter_t& time_usage){
            auto comp_start = timer::now();
            sdsl::load_from_file(m_ms, ms_path);
            time_usage.register_now("load_ms", comp_start);

            auto ds_start = timer::now();
            m_ms_sel = ms_sel_1_type(&m_ms);
            time_usage.register_now("select_init", ds_start);
        }

        void _show_vec(){
            cout << endl;
            for(int i=0; i<m_ms.size(); i++)
                cout << m_ms[i] << " ";
            cout << endl;
        }

        /* walk all the bits from bit_from to bit_to */
        size_type trivial_range_max(const size_type int_from, const size_type int_to) const {
            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, max_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;

            if(int_from > 0){
                bit_from = m_ms_sel(int_from);
                prev_ms = bit_from - 2 * (int_from - 1);
                i = bit_from + 1;
            }
            while (cnt1 < (int_to - int_from)) {
                if (m_ms[i] == 1) {
                    cur_ms = prev_ms + cnt0 - 1;
                    max_ms = std::max(max_ms, cur_ms);
                    prev_ms = cur_ms;
                    cnt0 = 0;
                    cnt1 += 1;
                } else {
                    cnt0 += 1;
                }
                i += 1;
            }
            return max_ms;
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
                    cum_ms = std::max(cum_ms, ms_value);
                    one_cnt += 1;
                }
                if (ms_idx and (ms_idx + 1) % block_size == 0) {
                    out_vec[out_idx++] = cum_ms;
                    cum_ms = 0;
                }
            }
        }
    };

    template<typename size_type>
    class none_partial_max_vector : sdsl_partial_max_vector<sdsl::bit_vector, sdsl::bit_vector::select_1_type, size_type> {
    public:
        typedef sdsl_partial_max_vector<sdsl::bit_vector, sdsl::bit_vector::select_1_type, size_type>  base_cls;
        typedef typename base_cls::counter_t counter_t;

        none_partial_max_vector(const string& ms_path) : base_cls(ms_path) {}

        none_partial_max_vector(const string& ms_path, counter_t& time_usage) : base_cls(ms_path, time_usage) {}

        size_type check_range_max(const size_type int_from, const size_type int_to) const {
            return base_cls::trivial_range_max(int_from, int_to);
        }

        size_type noindex_range_max(const size_type  int_from, const size_type int_to, const IndexedAlgorithm algo) const {
            if (int_from >= int_to)
                return 0;

            if(algo == IndexedAlgorithm::djamal){
                size_type bit_from = this->m_ms_sel(int_from + 1);
                size_type bit_to = this->m_ms_sel(int_to);
                size_type result_idx = 0;
                return (size_type) ms_range_max_fast64(int_from, bit_from, bit_to, this->m_ms.data(), &result_idx);
            } else if (algo == IndexedAlgorithm::trivial) {
                return base_cls::trivial_range_max(int_from, int_to);
            }
        }

        size_type extreme_one(size_type from, size_type to, bool right) const{
            if(right){
                while(this->m_ms[--to] == 0 and to >= from);
            } else {
                while(this->m_ms[from++] != 0 and from < to);
            }
            size_type idx = (right ? to : from);
            assert(this->m_ms[idx]);
            return idx;
        }

        size_type trivial_right_max(size_type ms_i, size_type bit_from, size_type bit_to) const {
            size_type max = 0;
            for(size_type i = bit_from; i <= bit_to; i++){
                if(this->m_ms[i]){
                    max = std::max(max, i - 2 * ms_i);
                    ms_i += 1;
                }
            }
            return max;
        }
        size_type trivial_left_max(size_type ms_i, size_type bit_from, size_type bit_to) const {
            assert(bit_from <= bit_to);
            size_type max = 0;
            for(size_type i = bit_to; i >= bit_from; i--){
                if(this->m_ms[i]){
                    max = std::max(max, i - 2 * ms_i);
                    ms_i -= 1;
                }
            }
            return max;
        }

        size_type indexed_range_max(const sdsl::int_vector<64>& ridx, const size_type int_from, const size_type int_to, const size_type bsize,
                const IndexedAlgorithm algo) const {
            if(algo == IndexedAlgorithm::djamal)
                throw string{"Not supported"};

            size_type _max = 0;
            size_type bit_from = this->m_ms_sel(int_from + 1);
            size_type bit_to = this->m_ms_sel(int_to);
            // k = bit_from / bsize
            // aligned_bit_from = select(rank(k + 1) + 1)
            size_type block_from = (bit_from / bsize);
            size_type block_from_inside = block_from + (bit_from % bsize > 0);
            size_type block_to = ((bit_to + 0) / bsize);
            size_type block_to_inside = block_to - ((bit_to + 1) % bsize > 0);
            //TODO: put a RMQ here
            //sdsl::range_maximum_sct<> ridx_rmq(&ridx);
            //size_type _block_max = ridx[ridx_rmq(block_from_inside, block_to_inside)];
            for(int i = block_from_inside; i <= block_to_inside; i++)
                _max = std::max(_max, ridx[i]);

            if(block_from < block_from_inside){
                assert(block_from + 1 == block_from_inside);
                { /*
                   * find max in block_from[bit_from ... k]
                   * if a 0-run runs over the block, k = index of first 1 in block_from_inside
                   * else k = idx of last 1 in block_from
                   */
                    size_type ms_i = int_from;
                    for(size_type i = bit_from; i < block_from_inside * bsize; i++){
                        if(this->m_ms[i]){
                            _max = std::max(_max, i - 2 * ms_i);
                            ms_i += 1;
                        }
                    }
                }
            }
            if(block_to > block_to_inside){
                assert(block_to == block_to_inside + 1);
                { /*
                   * find max in ms[k..bit_to]
                   * if a 0-run runs over the block, k = index of last 1 in block_to_inside + 1
                   * else k = idx of first 0 in block_to
                   */
                    size_type i = bit_to, ms_i = int_to - 1;
                    while(i >= bit_from){
                        if(this->m_ms[i]){
                            if(i / bsize < block_to)
                                break;
                            _max = std::max(_max, i - 2 * ms_i);
                            ms_i -= 1;
                        }
                        i -= 1;
                    }
                }
            }
            return _max;
        }
    };

    template <typename size_type>
    class rrr_partial_max_vector : sdsl_partial_max_vector<sdsl::rrr_vector<>, sdsl::rrr_vector<>::select_1_type, size_type>{
        inline uint64_t uncompress64(const size_type word_from) const {
            uint64_t bit_from = word_from * 64;
            uint8_t width = std::min<uint64_t>(64, this->m_ms.size() - bit_from);
            return this->m_ms.get_int(bit_from, width);
        }
    public:
        typedef sdsl_partial_max_vector<sdsl::rrr_vector<>, sdsl::rrr_vector<>::select_1_type, size_type>  base_cls;
        typedef typename base_cls::counter_t counter_t;

        rrr_partial_max_vector(const string& ms_path) : base_cls(ms_path) {}

        rrr_partial_max_vector(const string& ms_path, counter_t& time_usage) : base_cls(ms_path, time_usage) {}

        size_type check_range_max(const size_type int_from, const size_type int_to) const {
            return base_cls::trivial_range_max(int_from, int_to);
        }

        size_type noindex_range_max(const size_type  int_from, const size_type int_to, const IndexedAlgorithm algo) const {
            if (int_from >= int_to)
                return 0;

            if(algo == IndexedAlgorithm::djamal){
                size_type bit_from = this->m_ms_sel(int_from + 1);
                size_type bit_to = this->m_ms_sel(int_to);
                size_type prev_ms = bit_from - 2 * int_from;
                throw string{"not supported yet"};
                //return _bit_djamal_range_sum_fast(bit_from, bit_to, prev_ms);
            } else if (algo == IndexedAlgorithm::trivial) {
                return base_cls::trivial_range_max(int_from, int_to);
            }
        }
    };
}
#endif /* PARTIAL_MAX_VECTOR_HPP */
