#ifndef PARTIAL_MAX_VECTOR_HPP
#define PARTIAL_MAX_VECTOR_HPP

#include <bitset>
#include <string>
#include <iostream>
#include <fstream>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>
extern "C" {
    #include "virtual_smsb/ms_range_max.h"
    #include "virtual_smsb/naive_ms_range_max.h"
}
#include "counter.hpp"
#include "range_query.hpp"
#include "partial_op_vector.hpp"


namespace fdms {
    /* rle based class */
    template<typename vec_type, typename it_type, typename size_type>
    class rle_partial_max_vector : public rle_partial_op_vector<vec_type, it_type, size_type> {
    public:
        typedef rle_partial_op_vector<vec_type, it_type, size_type>  base_cls;
        typedef typename base_cls::rqres_t rqres_t;

        rle_partial_max_vector(const vec_type& v, it_type* it) : base_cls(v, it) {}

        rqres_t trivial(const size_type int_from, const size_type int_to) {
            base_cls::check_range(int_from, int_to);

            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, max_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;
            size_type max_ms_idx = int_from;

            if(int_from > 0){
                bit_from = base_cls::_select_0based(int_from - 1); //base_cls::m_it->select(int_from - 1);
                prev_ms = base_cls::_int_ms_0based(bit_from, int_from - 1); //bit_from - 2 * (int_from - 1);
                i = bit_from + 1;
            }
            while (cnt1 < (int_to - int_from)) {
                if (base_cls::m_it->isSet(i)) {
                    cur_ms = prev_ms + cnt0 - 1;
                    if(cur_ms > max_ms){
                        max_ms = cur_ms;
                        max_ms_idx = int_from + cnt1;
                    }
                    prev_ms = cur_ms;
                    cnt0 = 0;
                    cnt1 += 1;
                } else {
                    cnt0 += 1;
                }
                i += 1;
            }
            return rqres_t(max_ms_idx, max_ms);
        }

        rqres_t __djamal_fast(const size_type n_ones,
                const size_type bit_from, size_type prev_ms){
            size_type max_ms = 0, cnt1 = 0, i = bit_from;
            size_type max_ms_idx = 0; // means that the index is relative

            while (cnt1 < n_ones) {
                size_type j = base_cls::m_it->selectNext();
                size_type cnt0 = (j - i - 1);
                size_type cur_ms = prev_ms + cnt0 - 1;
                if(cur_ms > max_ms){
                    max_ms = cur_ms;
                    max_ms_idx = cnt1;
                }
                prev_ms = cur_ms;
                cnt1 += 1;
                i = j;
            }
            return rqres_t(max_ms_idx, max_ms);
        }

        rqres_t __djamal_faster(const size_type n_ones,
                const size_type bit_from, size_type prev_ms, const size_type ms_size){
            size_type max_ms = 0, cnt1 = 0, i = bit_from;
            size_type max_ms_idx = 0; // means that the index is relative
            std::pair<size_type, size_type> run_state;

            while(cnt1 < n_ones) {
                // ms[run_state.1] = 1 && there are run_state.2 1s after that
                run_state = base_cls::m_it->selectNextRun(ms_size);
                size_type cnt0 = run_state.first - i - 1;
                size_type cur_ms = prev_ms + cnt0 - 1;
                size_t limit = std::min(run_state.second + 1, n_ones - cnt1);

                for(size_type j = 0; j < limit; j++){
                    cur_ms = prev_ms + cnt0 - 1;
                    if(cur_ms > max_ms){
                        max_ms = cur_ms;
                        max_ms_idx = cnt1;
                    }
                    prev_ms = cur_ms;
                    cnt1 += 1;

                    if(cnt1 >= n_ones)
                        break;
                    cnt0 = 0;
                }
                i = run_state.first + run_state.second;
            }
            return rqres_t(max_ms_idx, max_ms);
        }

        rqres_t __djamal_fastest(const size_type n_ones,
                const size_type bit_from, size_type prev_ms, const size_type ms_size){
            size_type max_ms = 0, cnt1 = 0, i = bit_from;
            size_type max_ms_idx = 0; // means that the index is relative
            std::pair<size_type, size_type> run_state;

            while(cnt1 < n_ones) {
                // ms[run_state.1] = 1 && there are run_state.2 1s after that
                run_state = base_cls::m_it->selectNextRun(ms_size);
                size_type cnt0 = run_state.first - i - 1;
                size_type cur_ms = prev_ms + cnt0 - 1;
                size_t limit = std::min(run_state.second + 1, n_ones - cnt1);

                // the MS corresponding to the first 1 is max
                if(cur_ms > max_ms){
                    max_ms = cur_ms;
                    max_ms_idx = cnt1;
                }

                // the rest of 1s in this run -> smaller MS values
                cnt1 += limit;
                prev_ms = cur_ms - limit + 1;
                i = run_state.first + run_state.second;
            }
            return rqres_t(max_ms_idx, max_ms);
        }

        rqres_t djamal(const size_type int_from, const size_type int_to) {
            base_cls::check_range(int_from, int_to);

            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;
            rqres_t max_ms;

            if(int_from > 0){
                bit_from = base_cls::_select_0based(int_from - 1); //base_cls::m_it->select(int_from - 1);
                prev_ms = base_cls::_int_ms_0based(bit_from, int_from - 1); //bit_from - 2 * (int_from - 1);
                i = bit_from;
            } else {
                base_cls::m_it->select(0); // this will initialize the iterator
            }

            //max_ms = __djamal_fast(int_to - int_from, bit_from, prev_ms);
            //max_ms = __djamal_faster(int_to - int_from, bit_from, prev_ms, base_cls::m_ms.getSize());
            max_ms = __djamal_fastest(int_to - int_from, bit_from, prev_ms, base_cls::m_ms.getSize());
            return rqres_t(int_from + max_ms.index, max_ms.value);
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
            throw string("Indexed max on rle vectors not implemented");
        }
    };

    /* sdsl based class */
    template<typename vec_type, typename ms_sel_1_type, typename ms_rank_1_type, typename size_type>
    class sdsl_partial_max_vector : public sdsl_partial_op_vector<vec_type, ms_sel_1_type, ms_rank_1_type, sdsl::rmq_succinct_sct<false>, size_type>{
    protected:
        typedef typename sdsl::rmq_succinct_sct<false> idx_vector_t;
        typedef sdsl_partial_op_vector<vec_type, ms_sel_1_type, ms_rank_1_type, idx_vector_t, size_type>  base_cls;

    public:
        typedef typename base_cls::counter_t counter_t;
        typedef typename base_cls::rqres_t rqres_t;
        typedef pair<rqres_t, size_type> pair_t;

        ms_rank_1_type m_rb;

        sdsl_partial_max_vector(const string& ms_path) : base_cls(ms_path) {
            m_rb = ms_rank_1_type(&(base_cls::m_ms));

            for(auto k: {"range.int", "range.bit"})
                 base_cls::m_time_usage.reg[k] = static_cast<size_type>(0);
        }

        pair_t bit_trivial_shift(size_type int_from, size_type bit_idx, size_type bit_idx_stop) const {
            assert(base_cls::m_ms[bit_idx]);
            size_type _cur = 0, _max = 0, cnt1 = 0, _max_index = 0, ms_idx = int_from;
            do{
                if(base_cls::m_ms[bit_idx]){
                    _cur = bit_idx - 2 * ms_idx;
                    if(_cur > _max){
                        _max = _cur;
                        _max_index = int_from + cnt1;
                    }
                    ms_idx += 1;
                    cnt1 += 1;
                }
            } while(bit_idx++ < bit_idx_stop);
            return std::make_pair(rqres_t(_max_index, _max), cnt1);
        }

        rqres_t trivial(const size_type int_from, const size_type int_to) const {
            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, max_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;
            size_type max_ms_idx = int_from;

            if(int_from > 0){
                bit_from = base_cls::m_ms_sel(int_from);
                prev_ms = bit_from - 2 * (int_from - 1);
                i = bit_from + 1;
            }
            while (cnt1 < (int_to - int_from)) {
                if (base_cls::m_ms[i] == 1) {
                    cur_ms = prev_ms + cnt0 - 1;
                    if(cur_ms > max_ms){
                        max_ms = cur_ms;
                        max_ms_idx = int_from + cnt1;
                    }
                    prev_ms = cur_ms;
                    cnt0 = 0;
                    cnt1 += 1;
                } else {
                    cnt0 += 1;
                }
                i += 1;
            }
            return rqres_t(max_ms_idx, max_ms);
        }

        static void dump(const string& ms_path, const size_type block_size) {
            sdsl::int_vector_buffer<1> ms(ms_path, std::ios::in);
            sdsl::int_vector_buffer<64> out_vec(InputSpec::rdix_fname(ms_path, block_size), std::ios::out);

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
    class none_partial_max_vector : public sdsl_partial_max_vector<sdsl::bit_vector, sdsl::bit_vector::select_1_type, sdsl::bit_vector::rank_1_type, size_type> {
    public:
        typedef sdsl_partial_max_vector<sdsl::bit_vector, sdsl::bit_vector::select_1_type, sdsl::bit_vector::rank_1_type, size_type>  base_cls;
        typedef typename base_cls::counter_t counter_t;
        typedef typename base_cls::rqres_t rqres_t;
        typedef typename base_cls::idx_vector_t idx_vector_t;
        typedef typename base_cls::pair_t pair_t;
    private:
        rqres_t max_in_block(const size_type block_idx, const size_type bsize) const {
            size_type first_one_idx = block_idx * bsize;
            while(this->m_ms[first_one_idx] == 0)
                first_one_idx += 1;
            assert(first_one_idx < (block_idx + 1) * bsize);
            size_type ms_i = base_cls::m_rb(first_one_idx);
            return base_cls::bit_trivial_shift(ms_i, first_one_idx, (block_idx + 1) * bsize - 1).first;
        }

    public:
        none_partial_max_vector(const string& ms_path) : base_cls(ms_path) {
            std::vector<string> _keys = {
                    "algorithm.trivial_case",
                    "algorithm.rmq_scan", "algorithm.rmq_query",
                    "algorithm.trivial_scan",
                    "algorithm.trivial_scan.1", "algorithm.trivial_scan.2"
            };
            for(auto k: _keys)
                base_cls::m_time_usage.reg[k] = static_cast<size_type>(0);
            for(auto k: {"range.int", "range.bit", "range.block", "range.i_block"})
                base_cls::m_time_usage.reg[k] = static_cast<size_type>(0);
        }


        rqres_t noindex(const size_type int_from, const size_type int_to, const RangeAlgorithm algo) {
            base_cls::check_range(int_from, int_to);

            if(algo == RangeAlgorithm::djamal){
                size_type bit_from = base_cls::_select_0based(int_from);
                size_type bit_to = base_cls::_select_0based(int_to - 1); // since it is a half-open interval
                size_type result_idx = 0;
                size_type _max = (size_type) ms_range_max_fast64(int_from, bit_from, bit_to, this->m_ms.data(), &result_idx);
                return rqres_t(result_idx, _max);
            } else if (algo == RangeAlgorithm::trivial) {
                return base_cls::trivial(int_from, int_to);
            }
            throw string{"Bad algorithm in non-indexed max range."};
        }

        rqres_t indexed(const idx_vector_t& rmq,
                const size_type int_from, const size_type int_to, const size_type bsize,
                const RangeAlgorithm algo) {
            base_cls::check_range(int_from, int_to);
            if(algo == RangeAlgorithm::djamal)
                throw string{"Djamal algorithm not supported, for this indexed max-range."};

            auto cs1 = timer::now();
            size_type _max = 0, _max_index = 0;
            size_type bit_from = this->m_ms_sel(int_from + 1);
            size_type bit_to = this->m_ms_sel(int_to);

            size_type block_from = (bit_from / bsize);
            size_type block_from_inside = block_from + (bit_from % bsize > 0 && (block_from + 1) * bsize < bit_to);
            size_type block_to = ((bit_to + 0) / bsize);
            size_type block_to_inside = block_to - ((bit_to + 1) % bsize > 0 && block_to * bsize > bit_from);
            base_cls::m_time_usage.register_add("algorithm.init", cs1);

            // get some stats
            base_cls::m_time_usage.reg["range.int"] += static_cast<size_type>(int_to - int_from);
            base_cls::m_time_usage.reg["range.bit"] += static_cast<size_type>(bit_to - bit_from);
            base_cls::m_time_usage.reg["range.block"] += static_cast<size_type>(block_to - block_from);
            if(block_from_inside < block_to_inside)
                base_cls::m_time_usage.reg["range.i_block"] += static_cast<size_type>(block_to_inside - block_from_inside);

            if(block_from_inside >= block_to_inside){ // 1 or less proper inside blocks
                auto cs3 = timer::now();
                rqres_t res = base_cls::trivial(int_from, int_to);
                base_cls::m_time_usage.register_add("algorithm.trivial_case", cs3);
                return res;
            }

            size_type nr1_inside_blocks = int_to - int_from;

            // there are more than 1 proper inside blocks
            auto cs4 = timer::now();
            auto cs41 = timer::now();
            // section up to the first inside block
            if(block_from < block_from_inside){
                assert(block_from + 1 == block_from_inside);
                assert(this->m_ms[bit_from] ==  1);
                pair_t res = base_cls::bit_trivial_shift(int_from, bit_from, block_from_inside * bsize);
                _max = res.first.value;
                _max_index = res.first.index;
                nr1_inside_blocks -= res.second;
                base_cls::m_time_usage.reg["algorithm.i.1"] += static_cast<size_type>(block_from_inside * bsize - bit_from + 1);
            }
            base_cls::m_time_usage.register_add("algorithm.trivial_scan.1", cs41);
            // section following the last inside block
            if(block_to > block_to_inside){
                auto cs42 = timer::now();
                assert(block_to == block_to_inside + 1);
                size_type i = bit_to, ms_i = int_to - 1, block_limit = block_to * bsize;
                while(i >= bit_from and i >= block_limit){
                    if(this->m_ms[i]){
                        size_type _cur = i - 2 * ms_i;
                        if(_cur > _max){
                            _max = _cur;
                            _max_index = ms_i;
                        }
                        ms_i -= 1;
                        assert(nr1_inside_blocks > 0);
                        nr1_inside_blocks -= 1;
                    }
                    i -= 1;
                }
                base_cls::m_time_usage.register_add("algorithm.trivial_scan.2", cs42);
            }
            base_cls::m_time_usage.register_add("algorithm.trivial_scan", cs4);

            if(nr1_inside_blocks == 0)
                return rqres_t(0, _max);

            // there are some 1s in the proper inside blocks
            auto cs5 = timer::now();
            size_type block_idx = rmq(block_from_inside, block_to_inside);
            assert(block_from_inside <= block_idx and block_idx <= block_to_inside);
            base_cls::m_time_usage.register_add("algorithm.rmq_query", cs5);

            auto cs6 = timer::now();
            rqres_t proper_res = max_in_block(block_idx, bsize);
            rqres_t inproper_res = rqres_t(_max_index, _max);
            base_cls::m_time_usage.register_add("algorithm.rmq_scan", cs6);
            if(proper_res.value > inproper_res.value)
                return proper_res;
            else if (proper_res.value < inproper_res.value)
                return inproper_res;
            return rqres_t(std::min(proper_res.index, inproper_res.index), proper_res.value);
        }
    };

    template <typename size_type>
    class rrr_partial_max_vector : public sdsl_partial_max_vector<sdsl::rrr_vector<>, sdsl::rrr_vector<>::select_1_type, sdsl::rrr_vector<>::rank_1_type, size_type>{
    public:
        typedef sdsl_partial_max_vector<sdsl::rrr_vector<>, sdsl::rrr_vector<>::select_1_type, sdsl::rrr_vector<>::rank_1_type, size_type>  base_cls;
        typedef typename base_cls::counter_t counter_t;
        typedef typename base_cls::rqres_t rqres_t;
        typedef typename base_cls::idx_vector_t idx_vector_t;
        typedef typename base_cls::pair_t pair_t;
    private:
        inline uint64_t uncompress64(const size_type word_from) const {
            uint64_t bit_from = word_from * 64;
            uint8_t width = std::min<uint64_t>(64, this->m_ms.size() - bit_from);
            return this->m_ms.get_int(bit_from, width);
        }
        rqres_t max_in_block(const size_type block_idx, const size_type bsize) const {
            size_type first_one_idx = block_idx * bsize;
            while(this->m_ms[first_one_idx] == 0)
                first_one_idx += 1;
            assert(first_one_idx < (block_idx + 1) * bsize);
            size_type ms_i = base_cls::m_rb(first_one_idx);
            return base_cls::bit_trivial_shift(ms_i, first_one_idx, (block_idx + 1) * bsize - 1).first;
        }

    public:
        rrr_partial_max_vector(const string& ms_path) : base_cls(ms_path) {}

        rqres_t noindex(const size_type  int_from, const size_type int_to, const RangeAlgorithm algo) {
            base_cls::check_range(int_from, int_to);

            if(algo == RangeAlgorithm::djamal)
                throw string{"Not supported on non-indexed range-max rrr-vector."};

            return base_cls::trivial(int_from, int_to);
        }

        rqres_t indexed(const idx_vector_t& rmq,
                const size_type int_from, const size_type int_to, const size_type bsize,
                const RangeAlgorithm algo) {
            base_cls::check_range(int_from, int_to);
            if(algo == RangeAlgorithm::djamal)
                throw string{"Not supported on indexed range-max rrr-vector."};

            auto cs1 = timer::now();
            size_type _max = 0, _max_index = 0;
            size_type bit_from = this->m_ms_sel(int_from + 1);
            size_type bit_to = this->m_ms_sel(int_to);

            size_type block_from = (bit_from / bsize);
            size_type block_from_inside = block_from + (bit_from % bsize > 0 && (block_from + 1) * bsize < bit_to);
            size_type block_to = ((bit_to + 0) / bsize);
            size_type block_to_inside = block_to - ((bit_to + 1) % bsize > 0 && block_to * bsize > bit_from);
            base_cls::m_time_usage.register_add("algorithm.init", cs1);

            // get some stats
            base_cls::m_time_usage.reg["range.int"] += static_cast<size_type>(int_to - int_from);
            base_cls::m_time_usage.reg["range.bit"] += static_cast<size_type>(bit_to - bit_from);
            base_cls::m_time_usage.reg["range.block"] += static_cast<size_type>(block_to - block_from);
            if(block_from_inside < block_to_inside)
                base_cls::m_time_usage.reg["range.i_block"] += static_cast<size_type>(block_to_inside - block_from_inside);

            if(block_from_inside >= block_to_inside){ // 1 or less proper inside blocks
                auto cs3 = timer::now();
                rqres_t res = base_cls::trivial(int_from, int_to);
                base_cls::m_time_usage.register_add("algorithm.trivial_case", cs3);
                return res;
            }

            size_type nr1_inside_blocks = int_to - int_from;

            // there are more than 1 proper inside blocks
            auto cs4 = timer::now();
            auto cs41 = timer::now();
            // section up to the first inside block
            if(block_from < block_from_inside){
                assert(block_from + 1 == block_from_inside);
                assert(this->m_ms[bit_from] ==  1);
                pair_t res = base_cls::bit_trivial_shift(int_from, bit_from, block_from_inside * bsize);
                _max = res.first.value;
                _max_index = res.first.index;
                nr1_inside_blocks -= res.second;
                base_cls::m_time_usage.reg["algorithm.i.1"] += static_cast<size_type>(block_from_inside * bsize - bit_from + 1);
            }
            base_cls::m_time_usage.register_add("algorithm.trivial_scan.1", cs41);
            // section following the last inside block
            if(block_to > block_to_inside){
                auto cs42 = timer::now();
                assert(block_to == block_to_inside + 1);
                size_type i = bit_to, ms_i = int_to - 1, block_limit = block_to * bsize;
                while(i >= bit_from and i >= block_limit){
                    if(this->m_ms[i]){
                        size_type _cur = i - 2 * ms_i;
                        if(_cur > _max){
                            _max = _cur;
                            _max_index = ms_i;
                        }
                        ms_i -= 1;
                        assert(nr1_inside_blocks > 0);
                        nr1_inside_blocks -= 1;
                    }
                    i -= 1;
                }
                base_cls::m_time_usage.register_add("algorithm.trivial_scan.2", cs42);
            }
            base_cls::m_time_usage.register_add("algorithm.trivial_scan", cs4);

            if(nr1_inside_blocks == 0)
                return rqres_t(_max_index, _max);

            // there are some 1s in the inside blocks
            auto cs5 = timer::now();
            size_type block_idx = rmq(block_from_inside, block_to_inside);
            assert(block_from_inside <= block_idx and block_idx <= block_to_inside);
            base_cls::m_time_usage.register_add("algorithm.rmq_query", cs5);

            auto cs6 = timer::now();
            rqres_t proper_res = max_in_block(block_idx, bsize);
            rqres_t inproper_res = rqres_t(_max_index, _max);
            base_cls::m_time_usage.register_add("algorithm.rmq_scan", cs6);
            if(proper_res.value > inproper_res.value)
                return proper_res;
            else if (proper_res.value < inproper_res.value)
                return inproper_res;
            return rqres_t(std::min(proper_res.index, inproper_res.index), proper_res.value);
        }
    };
}
#endif /* PARTIAL_MAX_VECTOR_HPP */
