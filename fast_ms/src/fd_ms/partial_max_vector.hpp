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

namespace fdms {
    /* rle based class */
    template<typename vec_type, typename it_type, typename size_type>
    class rle_partial_max_vector {
    public:

        const vec_type& m_ms;
        it_type *m_it;

        rle_partial_max_vector(const vec_type& v, it_type* it) : m_ms{v}, m_it{it} {}

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

        size_type trivial(const size_type int_from, const size_type int_to) {
            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, max_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;

            if(int_from > 0){
                bit_from = m_it->select(int_from - 1);
                //cout << "+ " << int_from << " -> " << bit_from << endl;
                prev_ms = bit_from - 2 * (int_from - 1);
                i = bit_from + 1;
            }
            while (cnt1 < (int_to - int_from)) {
                if (m_it->isSet(i)) {
                    //(cerr << "MS[" << cnt1 - 1 << "] = " << prev_ms << ", SUM = " << max_ms << endl);
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

        size_type __djamal_fast(const size_type n_ones, const size_type bit_from, size_type prev_ms){
            size_type max_ms = 0, cnt1 = 0, i = bit_from;

            while (cnt1 < n_ones) {
                size_type j = m_it->selectNext();
                size_type cnt0 = (j - i - 1);
                size_type cur_ms = prev_ms + cnt0 - 1;
                max_ms = std::max(max_ms, cur_ms);
                prev_ms = cur_ms;
                cnt1 += 1;
                i = j;
            }
            return max_ms;
        }

        size_type __djamal_faster(const size_type n_ones,
                                         const size_type bit_from, size_type prev_ms, const size_type ms_size){
            size_type max_ms = 0, cnt1 = 0, i = bit_from;
            std::pair<size_type, size_type> run_state;

            while(cnt1 < n_ones) {
                //ms[run_state.1] = 1 && there are run_state.2 1s after that
                run_state = m_it->selectNextRun(ms_size);
                size_type cnt0 = run_state.first - i - 1;
                size_type cur_ms = prev_ms + cnt0 - 1;
                size_t limit = std::min(run_state.second + 1, n_ones - cnt1);

                for(size_type j = 0; j < limit; j++){
                    cur_ms = prev_ms + cnt0 - 1;
                    max_ms = std::max(max_ms, cur_ms);
                    prev_ms = cur_ms;
                    cnt1 += 1;

                    if(cnt1 >= n_ones)
                        break;
                    cnt0 = 0;
                }
                i = run_state.first + run_state.second;
            }
            return max_ms;
        }

        size_type __djamal_fastest(const size_type n_ones, const size_type bit_from, size_type prev_ms, const size_type ms_size){
            size_type max_ms = 0, cnt1 = 0, i = bit_from;
            std::pair<size_type, size_type> run_state;

            while(cnt1 < n_ones) {
                //ms[run_state.1] = 1 && there are run_state.2 1s after that
                run_state = m_it->selectNextRun(ms_size);
                size_type cnt0 = run_state.first - i - 1;
                size_type cur_ms = prev_ms + cnt0 - 1;
                size_t limit = std::min(run_state.second + 1, n_ones - cnt1);

                // the MS corresponding to the first 1 is max
                max_ms = std::max(max_ms, cur_ms);

                // the rest of 1s in this run -> smaller MS values
                cnt1 += limit;
                prev_ms = cur_ms - limit + 1;
                i = run_state.first + run_state.second;
            }
            return max_ms;
        }

        size_type djamal(const size_type int_from, const size_type int_to){
            size_type bit_from = 0;
            size_type prev_ms = 1, cur_ms = 0, max_ms = 0;
            size_type cnt1 = 0, cnt0 = 0, i = bit_from;

            if(int_from > 0){
                bit_from = m_it->select(int_from - 1);
                //cout << "+ " << int_from << " -> " << bit_from << endl;
                prev_ms = bit_from - 2 * (int_from - 1);
                i = bit_from;
            } else {
                m_it->select(0); // this will initialize the iterator
            }

            //max_ms = __djamal_fast(int_to - int_from, bit_from, prev_ms);
            //max_ms = __djamal_faster(int_to - int_from, bit_from, prev_ms, m_ms.getSize());
            max_ms = __djamal_fastest(int_to - int_from, bit_from, prev_ms, m_ms.getSize());
            return max_ms;
        }

        size_type indexed(sdsl::int_vector<64>& ridx,
                const size_type from, const size_type to, const size_type bsize) {

            assert(from < to);
            throw string("Indexed max on rle vectors not implemented");
        }
    };

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
    class none_partial_max_vector : public sdsl_partial_max_vector<sdsl::bit_vector, sdsl::bit_vector::select_1_type, size_type> {
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

        size_type indexed_range_max(const sdsl::int_vector<64>& ridx, const sdsl::rmq_succinct_sct<false> &rmq,
                const sdsl::bit_vector::rank_1_type &rb,
                size_type int_from, const size_type int_to, const size_type bsize,
                const IndexedAlgorithm algo) const {
            if(algo == IndexedAlgorithm::djamal)
                throw string{"Not supported"};

            size_type _max = 0;
            size_type bit_from = this->m_ms_sel(int_from + 1);
            size_type bit_to = this->m_ms_sel(int_to);
            // k = bit_from / bsize
            // aligned_bit_from = select(rank(k + 1) + 1)
            size_type block_from = (bit_from / bsize);
            size_type block_from_inside = block_from + (bit_from % bsize > 0 && (block_from + 1) * bsize < bit_to);
            size_type block_to = ((bit_to + 0) / bsize);
            size_type block_to_inside = block_to - ((bit_to + 1) % bsize > 0 && block_to * bsize > bit_from);

            if(block_from_inside >= block_to_inside){ // 1 or less proper inside blocks
                return base_cls::trivial_range_max(int_from, int_to);
                //TODO: the method below is faster!!!
                size_type ms_i = int_from;
                for(size_type i = bit_from; i <= bit_to; i++){
                    if(this->m_ms[i]){
                        _max = std::max(_max, i - 2 * ms_i);
                        ms_i += 1;
                    }
                }
                return _max;
            }
            // there are 1 or more proper inside blocks
            size_type block_idx = rmq(block_from_inside, block_to_inside);
            assert(block_from_inside <= block_idx and block_idx <= block_to_inside);
            size_type first_one_idx = block_idx * bsize;
            {
                while(this->m_ms[first_one_idx] == 0)
                    first_one_idx += 1;
                assert(first_one_idx < (block_idx + 1) * bsize);
            }
            {
                size_type ms_i = rb(first_one_idx), i = first_one_idx;
                do{
                    if(this->m_ms[i]){
                        _max = std::max(_max, i - 2 * ms_i);
                        ms_i += 1;
                    }
                } while(++i < (block_idx + 1) * bsize);
            }

            if(block_from < block_from_inside){
                assert(block_from + 1 == block_from_inside);
                size_type ms_i = int_from;
                for(size_type i = bit_from; i < block_from_inside * bsize; i++){
                    if(this->m_ms[i]){
                        _max = std::max(_max, i - 2 * ms_i);
                        ms_i += 1;
                    }
                }
            }

            if(block_to > block_to_inside){
                assert(block_to == block_to_inside + 1);
                size_type i = bit_to, ms_i = int_to - 1;
                while(i >= bit_from){
                    if(this->m_ms[i]){
                        if(i / bsize < block_to)
                            break;
                        _max = std::max(_max, i - 2 * ms_i);
                        ms_i -= 1;
                    }
                    if(i == 0)
                        break;
                    i -= 1;
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
