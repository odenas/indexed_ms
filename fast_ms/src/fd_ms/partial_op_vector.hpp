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
