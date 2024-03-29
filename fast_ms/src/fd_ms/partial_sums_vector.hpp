#ifndef PARTIAL_SUMS_VECTOR_HPP
#define PARTIAL_SUMS_VECTOR_HPP

#include <bitset>
#include <string>
#include <iostream>
#include <fstream>

#include <sdsl/vectors.hpp>
extern "C" {
    #include "virtual_smsb/ms_range_sum.h"
    #include "virtual_smsb/naive_ms_range_sum.h"
}
#include "counter.hpp"
#include "range_query.hpp"
#include "partial_op_vector.hpp"


namespace fdms {
    /* rle - based class */
    template<typename vec_type, typename it_type, typename size_type>
    class rle_partial_sums_vector : public rle_partial_op_vector<vec_type, it_type, size_type> {
    private:
        /**
         * compute the sum of terms: base + (base - 1) + ... (base - n_terms - 1)
         */
        static size_type sum_consecutive_ms(const size_type base, const size_type n_terms) {
            if (n_terms > base)
                throw string{"base = " + to_string(base) + "must be > n_terms (" + to_string(n_terms) + ")"};
            if(n_terms == 0)
                throw string{"n_terms must be positive. got " + to_string(n_terms)};

            size_type a = base * (n_terms);
            size_type b = (n_terms * (n_terms - 1)) / 2;
            assert(a > b);

            return (a - b);
        }


    public:
        typedef rle_partial_op_vector<vec_type, it_type, size_type>  base_cls;
        typedef typename base_cls::rqres_t rqres_t;

        rle_partial_sums_vector(const vec_type& v, it_type* it) : base_cls(v, it) {}


        rqres_t trivial(const size_type int_from, const size_type int_to) {
            base_cls::check_range(int_from, int_to);

            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, sum_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;

            if(int_from > 0){
                bit_from = base_cls::_select_0based(int_from - 1); //base_cls::m_it->select(int_from - 1);
                prev_ms = base_cls::_int_ms_0based(bit_from, int_from - 1); //bit_from - 2 * (int_from - 1);
                i = bit_from + 1;
            }
            while (cnt1 < (int_to - int_from)) {
                if (base_cls::m_it->isSet(i)) {
                    cur_ms = prev_ms + cnt0 - 1; // since MS_i - MS_{i-1} + 1 = nzeros
                    sum_ms += cur_ms;
                    prev_ms = cur_ms;
                    cnt0 = 0;
                    cnt1 += 1;
                } else {
                    cnt0 += 1;
                }
                i += 1;
            }
            return rqres_t(0, sum_ms);
        }

        rqres_t __djamal_fast(const size_type n_ones,
                const size_type bit_from, size_type prev_ms){
            size_type sum_ms = 0, cnt1 = 0, i = bit_from;

            while (cnt1 < n_ones) {
                size_type j = base_cls::m_it->selectNext();
                size_type cnt0 = (j - i - 1);
                size_type cur_ms = prev_ms + cnt0 - 1;
                sum_ms += cur_ms;
                prev_ms = cur_ms;
                cnt1 += 1;
                i = j;
            }
            return rqres_t(0, sum_ms);
        }

        rqres_t __djamal_faster(const size_type n_ones,
                const size_type bit_from, size_type prev_ms, const size_type ms_size){
            size_type sum_ms = 0, cnt1 = 0, i = bit_from;
            std::pair<size_type, size_type> run_state;

            while(cnt1 < n_ones) {
                run_state = base_cls::m_it->selectNextRun(ms_size);
                size_type cnt0 = run_state.first - i - 1;
                size_type cur_ms = prev_ms + cnt0 - 1;
                size_t limit = std::min(run_state.second + 1, n_ones - cnt1);

                //size_type sum_check = sum_ms + sum_consecutive_ms(prev_ms + cnt0 - 1, limit);
                for(size_type j = 0; j < limit; j++){
                    cur_ms = prev_ms + cnt0 - 1;
                    sum_ms += cur_ms;
                    prev_ms = cur_ms;
                    cnt1 += 1;

                    if(cnt1 >= n_ones)
                        break;
                    cnt0 = 0;
                }
                i = run_state.first + run_state.second;
            }
            return rqres_t(0, sum_ms);
        }

        rqres_t __djamal_fastest(const size_type n_ones,
                const size_type bit_from, size_type prev_ms, const size_type ms_size) {
            size_type sum_ms = 0, cnt1 = 0, i = bit_from;
            std::pair<size_type, size_type> run_state;

            while(cnt1 < n_ones) {
                run_state = base_cls::m_it->selectNextRun(ms_size);  // run_state: (index of 1st 1, run length)
                size_type cnt0 = run_state.first - i - 1;
                size_type cur_ms = prev_ms + cnt0 - 1;
                size_type limit = std::min(run_state.second + 1, n_ones - cnt1);

                try{
                    sum_ms += sum_consecutive_ms(cur_ms, limit);
                } catch (string s){
                    throw s;
                }
                cnt1 += limit;
                prev_ms = cur_ms - limit + 1;
                i = run_state.first + run_state.second;
            }
            return rqres_t(0, sum_ms);
        }

        rqres_t djamal(const size_type int_from, const size_type int_to) {
            base_cls::check_range(int_from, int_to);

            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;
            rqres_t sum_ms;

            if(int_from > 0){
                bit_from = base_cls::_select_0based(int_from - 1); //base_cls::m_it->select(int_from - 1);
                prev_ms = base_cls::_int_ms_0based(bit_from, int_from - 1); //bit_from - 2 * (int_from - 1);
                i = bit_from;
            } else {
                base_cls::m_it->select(0); // this will initialize the iterator
            }

            //sum_ms = __djamalfast(int_to - int_from, bit_from, prev_ms);
            //sum_ms = __djamal_faster(int_to - int_from, bit_from, prev_ms, base_cls::m_ms.getSize());
            sum_ms = __djamal_fastest(int_to - int_from, bit_from, prev_ms, base_cls::m_ms.getSize());
            return sum_ms;
        }

        rqres_t noindex(const size_type int_from, const size_type int_to, const RangeAlgorithm algo) {
            base_cls::check_range(int_from, int_to);

            if(algo == RangeAlgorithm::djamal)
                return djamal(int_from, int_to);
            if (algo == RangeAlgorithm::trivial)
                return trivial(int_from, int_to);
            throw string{"Bad algorithm."};
        }

        rqres_t indexed(sdsl::int_vector<64>& ridx,
                const size_type from, const size_type to, const size_type bsize) {

            base_cls::check_range(from, to);

            rqres_t to_sum = indexed_prefix(ridx, to - 1, bsize);
            rqres_t from_sum = (from == 0 ? rqres_t() : indexed_prefix(ridx, from - 1, bsize));
            assert(from_sum.value <= to_sum.value);
            return rqres_t(0, to_sum.value - from_sum.value);
        }

        rqres_t indexed_prefix(sdsl::int_vector<64>& ridx, const size_type to_ms_idx, const size_type bsize) {
            //cerr << "[indexed (" << flags.block_size << ")] " << to_ms_idx << endl;

            // index of last term of sum
            size_type int_ms_idx = to_ms_idx;
            size_type bit_ms_idx = base_cls::m_it->select(int_ms_idx);
            size_type block_idx = bit_ms_idx / bsize;
            size_type sum_ms = 0; // to be subtracted from ridx[block_idx]
            {
                size_type prev_ms = bit_ms_idx - (2 * int_ms_idx); // needed for 1st term beyond the sum
                size_type nzeros = 0;

                // loop from bit_ms_idx + 1 to the end of the block
                for (size_type i = bit_ms_idx + 1; i < (block_idx + 1) * bsize; i++) {
                    if (base_cls::m_it->isSet(i) == 1) {
                        size_type cur_ms = (prev_ms + nzeros - 1); // since MS_i - MS_{i-1} + 1 = nzeros
                        sum_ms += cur_ms;
                        prev_ms = cur_ms;
                        nzeros = 0;
                    } else {
                        nzeros += 1;
                    }
                }
            }
            size_type answer = ridx[block_idx] - sum_ms;
            return rqres_t(0, answer);
        }
    };

    /* sdsl based class */
    template<typename vec_type, typename ms_sel_1_type, typename ms_rank_1_type, typename size_type>
    class sdsl_partial_sums_vector : public sdsl_partial_op_vector<vec_type, ms_sel_1_type, ms_rank_1_type, sdsl::int_vector<64>, size_type>{
    protected:
        typedef sdsl::int_vector<64> idx_vector_t;
        typedef pair<size_type, size_type> pair_t;

        /**
         * compute the sum up to the given position.
         * return the bit_index of the position and the sum
         */
        pair_t _slow_indexed_prefix(const idx_vector_t& ridx, const size_type int_to, const size_type bsize) {
            auto comp_start = timer::now();
            size_type int_ms_idx = int_to;
            size_type bit_ms_idx = base_cls::m_ms_sel(int_ms_idx + 1);
            size_type block_idx = bit_ms_idx / bsize;
            base_cls::m_time_usage.register_add("algorithm.p1", comp_start);

            comp_start = timer::now();
            size_type sum_ms = 0; // to be subtracted from ridx[block_idx]
            {
                size_type prev_ms = bit_ms_idx - (2 * int_ms_idx); // needed for 1st term beyond the sum
                size_type nzeros = 0;

                // loop from bit_ms_idx + 1 to the end_of_block
                size_type end_of_block = std::min<size_type>(base_cls::m_ms.size(), (block_idx + 1) * bsize);
                for (size_type i = bit_ms_idx + 1; i < end_of_block; i++) {
                    if (base_cls::m_ms[i] == 1) {
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
            base_cls::m_time_usage.register_add("algorithm.p2", comp_start);
            return std::make_pair(bit_ms_idx, answer);
        }

    public:
        typedef sdsl_partial_op_vector<vec_type, ms_sel_1_type, ms_rank_1_type, idx_vector_t, size_type>  base_cls;
        typedef typename base_cls::rqres_t rqres_t;
        typedef typename base_cls::counter_t counter_t;

        sdsl_partial_sums_vector(const string& ms_path) : base_cls{ms_path} {
            for(auto k: {"range.int", "range.bit"})
                 base_cls::m_time_usage.reg[k] = static_cast<size_type>(0);

            std::vector<string> _keys = {
                "algorithm.p1", "algorithm.p2",
            };
            for(auto k: _keys)
                base_cls::m_time_usage.register_now(k, timer::now());
        }

        /* walk all the bits from bit_from to bit_to */
        rqres_t trivial(const size_type int_from, const size_type int_to) const {
            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, sum_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;

            if(int_from > 0){
                bit_from = base_cls::m_ms_sel(int_from);
                prev_ms = bit_from - 2 * (int_from - 1);
                i = bit_from + 1;
            }
            while (cnt1 < (int_to - int_from)) {
                if (base_cls::m_ms[i] == 1) {
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
            return rqres_t(0, sum_ms);
        }

        static void dump(const string& ms_path, const size_type block_size) {
            sdsl::int_vector_buffer<1> ms(ms_path, std::ios::in);
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

    template<typename size_type>
    class none_partial_sums_vector : public sdsl_partial_sums_vector<sdsl::bit_vector, sdsl::bit_vector::select_1_type, sdsl::bit_vector::rank_1_type, size_type> {
        typedef sdsl_partial_sums_vector<sdsl::bit_vector, sdsl::bit_vector::select_1_type, sdsl::bit_vector::rank_1_type, size_type>  base_cls;
        typedef typename base_cls::counter_t counter_t;
        typedef typename base_cls::idx_vector_t idx_vector_t;
        typedef typename base_cls::rqres_t rqres_t;
        typedef typename base_cls::pair_t pair_t;

        pair_t _indexed_prefix(const idx_vector_t& ridx, const size_type int_to, const size_type bsize,
                const RangeAlgorithm algo) {

            if(algo == RangeAlgorithm::trivial) {
                return this->_slow_indexed_prefix(ridx, int_to - 1, bsize);
            }
            else if (algo == RangeAlgorithm::djamal) {
                auto comp_start = timer::now();
                size_type bit_to = this->m_ms_sel(int_to + 1);
                size_type prev_ms = bit_to - 2 * int_to;
                size_type block_to = bit_to / bsize;
                size_type bit_end_of_block_to = std::min<size_type>(this->m_ms.size(), (block_to + 1) * bsize) - 1;
                base_cls::m_time_usage.register_add("algorithm.p1", comp_start);

                comp_start = timer::now();
                size_type sum_ms = range_ms_sum_fast64(prev_ms, bit_to, bit_end_of_block_to, this->m_ms.data());
                size_type answer = ridx[block_to] - sum_ms;
                base_cls::m_time_usage.register_add("algorithm.p2", comp_start);
                return std::make_pair(bit_to, answer);
            }
            throw string{"Bad algo."};
        }

    public:

        none_partial_sums_vector(const string& ms_path) : base_cls(ms_path) {}

        rqres_t noindex(const size_type  int_from, const size_type int_to, const RangeAlgorithm algo) {
            base_cls::check_range(int_from, int_to);

            if(algo == RangeAlgorithm::djamal){
                size_type bit_from = base_cls::_select_0based(int_from);
                size_type bit_to = base_cls::_select_0based(int_to - 1);  // since it's a half open interval
                size_type prev_ms = base_cls::_int_ms_0based(bit_from, int_from);
                return rqres_t(0, (size_type) range_ms_sum_fast64(prev_ms, bit_from, bit_to, this->m_ms.data()));
            } else if (algo == RangeAlgorithm::trivial) {
                return base_cls::trivial(int_from, int_to);
            }
            throw string{"Bad algorithm in non-indexed sum range."};
        }

        rqres_t indexed(const idx_vector_t& ridx, const size_type from, const size_type to, const size_type bsize,
                const RangeAlgorithm algo) {
            base_cls::check_range(from, to);

            base_cls::m_time_usage.reg["range.int"] += static_cast<size_type>(to - from);

            pair_t res_to = _indexed_prefix(ridx, to, bsize, algo);
            pair_t res_from = (from == 0 ? std::make_pair(static_cast<size_type>(0), static_cast<size_type>(0)) : _indexed_prefix(ridx, from, bsize, algo));

            base_cls::m_time_usage.reg["range.bit"] += static_cast<size_type>(res_to.first - res_from.first);
            assert(res_from.second <= res_to.second);
            return rqres_t(0, res_to.second - res_from.second);
        }
    };

    template<typename size_type>
    class rrr_partial_sums_vector : public sdsl_partial_sums_vector<sdsl::rrr_vector<>, sdsl::rrr_vector<>::select_1_type, sdsl::rrr_vector<>::rank_1_type, size_type> {
        typedef sdsl_partial_sums_vector<sdsl::rrr_vector<>, sdsl::rrr_vector<>::select_1_type, sdsl::rrr_vector<>::rank_1_type, size_type>  base_cls;
        typedef typename base_cls::counter_t counter_t;
        typedef typename base_cls::idx_vector_t idx_vector_t;
        typedef typename base_cls::rqres_t rqres_t;
        typedef typename base_cls::pair_t pair_t;

        inline uint64_t uncompress64(const size_type word_from) const {
            uint64_t bit_from = word_from * 64;
            uint8_t width = std::min<uint64_t>(64, this->m_ms.size() - bit_from);
            return this->m_ms.get_int(bit_from, width);
        }

        size_type _bit_djamal(size_type bit_from, size_type bit_to, const size_type prev_ms) const {
            size_type word_from = (bit_from + 1) / 64;
            size_type word_to = bit_to / 64;
            size_type ms_sum = prev_ms;
            size_type curr_ms_sum;
            size_type virtual_ms = prev_ms;
            uint64_t vbuff = 0;

            if(bit_to <= bit_from)
                return prev_ms;
            bit_from = (bit_from + 1) % 64;
            bit_to = bit_to % 64;

            vbuff = uncompress64(word_from);
            if(word_from == word_to) {
                range_ms_sum_fast_ext64(virtual_ms, bit_from, bit_to,
                        &vbuff, &virtual_ms, &curr_ms_sum);
                return ms_sum + curr_ms_sum;
            }

            range_ms_sum_fast_ext64(virtual_ms, bit_from, 63,
                    &vbuff, &virtual_ms, &curr_ms_sum);
            ms_sum += curr_ms_sum;
            for(size_type i = word_from + 1; i < word_to; i++) {
                vbuff = uncompress64(i);
                range_ms_sum_fast_ext_word64(virtual_ms, vbuff,
                        &virtual_ms, &curr_ms_sum);
                ms_sum += curr_ms_sum;
            }
            vbuff = uncompress64(word_to);
            range_ms_sum_fast_ext64(virtual_ms, 0, bit_to,
                    &vbuff, &virtual_ms, &curr_ms_sum);
            return ms_sum + curr_ms_sum;
        }

        pair_t _indexed_prefix(const idx_vector_t& ridx, const size_type int_to, const size_type bsize,
                const RangeAlgorithm algo) {

            if(algo == RangeAlgorithm::trivial){
                return this->_slow_indexed_prefix(ridx, int_to -  1, bsize);
            } else if (algo == RangeAlgorithm::djamal) {
                size_type bit_to = this->m_ms_sel(int_to + 1);
                size_type prev_ms = bit_to - 2 * int_to;
                size_type block_to = bit_to / bsize;
                size_type sum_ms = _bit_djamal(bit_to, std::min<size_type>(this->m_ms.size(), (block_to + 1) * bsize) - 1, prev_ms);
                size_type partial_sum = ridx[block_to];
                return std::make_pair(bit_to, partial_sum - sum_ms);
            }
            throw string{"Bad algorithm in indexed sum range."};
        }

    public:
        rrr_partial_sums_vector(const string& ms_path) : base_cls(ms_path) {}

        /**
         * naive method that makes use of partial sums for queries [from_index, to_index)
         * by calling indexed_range_sum_prefix
         */
        rqres_t indexed(const idx_vector_t& ridx,  const size_type from, const size_type to, const size_type bsize,
                const RangeAlgorithm algo) {
            base_cls::check_range(from, to);

            base_cls::m_time_usage.reg["range.int"] += static_cast<size_type>(to - from);
            pair_t res_to = _indexed_prefix(ridx, to, bsize, algo);
            pair_t res_from = (from == 0 ? std::make_pair(static_cast<size_type>(0), static_cast<size_type>(0)) : _indexed_prefix(ridx, from, bsize, algo));

            base_cls::m_time_usage.reg["range.bit"] += static_cast<size_type>(res_to.first - res_from.first);
            assert(res_from.second <= res_to.second);
            return rqres_t(0, res_to.second - res_from.second);
        }

        rqres_t noindex(const size_type  int_from, const size_type int_to, const RangeAlgorithm algo) {
            base_cls::check_range(int_from, int_to);

            if(algo == RangeAlgorithm::djamal){
                size_type bit_from = base_cls::_select_0based(int_from);
                size_type bit_to = base_cls::_select_0based(int_to - 1); // since it's a half-open interval
                size_type prev_ms = base_cls::_int_ms_0based(bit_from, int_from);
                return rqres_t(0, _bit_djamal(bit_from, bit_to, prev_ms));
            } else if (algo == RangeAlgorithm::trivial) {
                return base_cls::trivial(int_from, int_to);
            }
            throw string{"Bad algorithm in non-indexed sum range."};
        }

    };
}
#endif /* PARTIAL_SUMS_VECTOR_HPP */
