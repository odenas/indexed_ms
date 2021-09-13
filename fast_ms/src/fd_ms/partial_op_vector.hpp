#ifndef PARTIAL_OP_VECTOR_HPP
#define PARTIAL_OP_VECTOR_HPP

#include <bitset>
#include <string>
#include <iostream>
#include <fstream>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>
#include "counter.hpp"
#include "range_query.hpp"

namespace fdms {

    /**
     * compute sel_i - 2i. check that i >= 0 and that the result is positive
     */
    template<typename size_type>
    size_type __int_ms_0based(const size_type sel_i, const size_type i) {
        assert(i >= 0);
        size_type j = 2 * i;
        if(sel_i < j){
            throw string{
                "i = " + to_string(i) +
                ", sel(i) = " + to_string(sel_i) +
                ", but 2i = " + to_string(j) +
                "\n 2i must be greater than sel(i)"
            };
        }
        return sel_i - j;
    }

    /* rle based class */
    template<typename vec_type, typename it_type, typename size_type>
    class rle_partial_op_vector {
    protected:
        typedef rq_result<size_type> rqres_t;

        const vec_type& m_ms;
        it_type *m_it;

        rle_partial_op_vector(const vec_type& v, it_type* it) : m_ms{v}, m_it{it} {}

        void check_range(const size_type from, const size_type to) const {
            if (from >= to){
                string range = string{"[" + to_string(from) + ", " + to_string(to) + ")"};
                throw string{"Empty range: " + range + "."};
            }
        }

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

        /**
        * compute sel_i - 2i. check that i >= 0 and that the result is positive
        */
        size_type _int_ms_0based(const size_type sel_i, const size_type i) const {
            try {
                return __int_ms_0based<size_type>(sel_i, i);
            } catch (string s){
                throw s;
            }
        }

        /**
        * compute MS[i+1] from MS[i] and n_zeros between the two using
        * MS_i - MS_{i-1} + 1 = nzeros between the two
        */
        size_type _next_int_ms(const size_type cur_ms, const size_type n_zeros) const {
            return cur_ms + n_zeros - 1;
        }

        size_type _select_0based(const size_type i) const {
            assert(i >= 0);
            size_type sel_i = m_it->select(i); // rle is 0-based
            if(!m_it->isSet(sel_i)){
                throw string{
                    "ms[select(" + to_string(i) + ") = " + to_string(sel_i) + "] = 0"
                };
            }
            return sel_i;
        }
    };


    /* sdsl based class */
    template<typename vec_type,
             typename ms_sel_1_type,
             typename ms_rank_1_type,
             typename idx_vector_t,
             typename size_type>
    class sdsl_partial_op_vector {
    protected:
        typedef sdsl::int_vector_buffer<1> buff_vec_t;
        typedef Counter<size_type> counter_t;
        typedef rq_result<size_type> rqres_t;

        size_type _int_ms_0based(size_type sel_i, size_type i) const {
            try {
                return __int_ms_0based<size_type>(sel_i, i);
            } catch (string s){
                throw s;
            }
        }

        size_type _select_0based(const size_type i) const {
            assert(i >= 0);
            size_type sel_i = m_ms_sel(i + 1); // sdsl is 1-based
            if(m_ms[sel_i] == 0){
                string select_message = {
                    "select(" + to_string(i) + ") = " + to_string(sel_i)
                };
                throw string{"ms[" + select_message + "] = 0"};
            }
            return sel_i;
        }


  public:
        vec_type m_ms;
        ms_sel_1_type m_ms_sel;
        counter_t m_time_usage;

        sdsl_partial_op_vector(const string& ms_path) {
            auto ds_start = timer::now();
            sdsl::load_from_file(m_ms, ms_path);
            m_ms_sel = ms_sel_1_type(&m_ms);
            m_time_usage.register_now("init", ds_start);
        }

        void check_range(const size_type from, const size_type to) const {
            if (from >= to)
                throw string{"Empty range: [" + std::to_string(from) + ", " + std::to_string(to) + ")."};
        }

        /**
         * compute MS[i] on the ms-bitvector as select(int_index) - 2*int_index
         */
        size_type int_ms(const size_type int_index){
            try{
                size_type bit_index = _select_0based(int_index);
                return _int_ms_0based(bit_index, int_index);
            } catch (string s){
                throw s;
            }
        }

        virtual rqres_t noindex(
            const size_type int_from, const size_type int_to,
            const RangeAlgorithm algo) = 0;

        virtual rqres_t indexed(
            const idx_vector_t& rmq,
            const size_type int_from, const size_type int_to, const size_type bsize,
            const RangeAlgorithm algo) = 0;
    };

}
#endif /* PARTIAL_OP_VECTOR_HPP */
