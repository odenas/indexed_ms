#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/select_support.hpp>

#include "fd_ms/counter.hpp"
#include "fd_ms/range_query.hpp"
#include "fd_ms/partial_op_vector.hpp"
#include "fd_ms/partial_sums_vector.hpp"
#include "fd_ms/partial_max_vector.hpp"

#include "rlcsa/bits/bitvector.h"
#include "rlcsa/bits/rlevector.h"
#include "rlcsa/bits/deltavector.h"
#include "rlcsa/bits/succinctvector.h"
#include "rlcsa/bits/nibblevector.h"

using namespace std;

namespace fdms {
    class rq_dispatcher{
    protected:
        typedef uint64_t size_type;
        static size_type random_index(const size_type max_idx) {
            return static_cast<size_type> (max_idx * static_cast<unsigned long> (std::rand()) / (RAND_MAX + 1UL));
        }
    public:
        typedef Counter<size_type> counter_t;
        typedef rq_result<size_type> rqres_t;
    };

    class sdsl_rq_dispatcher : rq_dispatcher {
        typedef pair<size_type, size_type> pair_t;

        static void _load_time_ridx(sdsl::int_vector<64>& ridx, const string ridx_path, counter_t& time_usage) {
            auto ds_start = timer::now();
            sdsl::load_from_file(ridx, ridx_path);
            time_usage.register_now("load_partial_sums", ds_start);
        }

        static void add_times_to(counter_t& to, counter_t& from){
            for(auto item: from.reg)
                to.reg[item.first] = item.second;
        }

    public:
        static rqres_t none_noindex(const string ms_path,
                const size_type from_idx, const size_type to_idx,
                const RangeAlgorithm algo, const RangeOperation op){
            counter_t  time_usage;
            if(op == RangeOperation::r_sum){
                return none_partial_sums_vector<size_type>(ms_path).noindex(from_idx, to_idx, algo);
            } else {
                return none_partial_max_vector<size_type>(ms_path).noindex(from_idx, to_idx, algo);
            }
        }
        template<typename op_cls>
        static void _noindex_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage, const RangeAlgorithm algo){

            op_cls psum(ms_path);
            auto comp_start = timer::now();
            for (int k = 0; k < nqueries; k++) {
                size_type start = random_index(from_idx_max);
                psum.noindex(start, start + range_size, algo);
                time_usage.register_now("algorithm", comp_start);
                add_times_to(time_usage, psum.m_time_usage);
            }
        }
        static void none_noindex_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage, const RangeAlgorithm algo, const RangeOperation op){

            if(op == RangeOperation::r_sum){
                _noindex_profile<none_partial_sums_vector<size_type>>(ms_path, nqueries, range_size, from_idx_max, time_usage, algo);
            } else {
                _noindex_profile<none_partial_max_vector<size_type>>(ms_path, nqueries, range_size, from_idx_max, time_usage, algo);
            }
        }

        static rqres_t rrr_noindex(const string ms_path,
                const size_type from_idx, const size_type to_idx,
                const RangeAlgorithm algo, const RangeOperation op){
            counter_t  time_usage;
            if(op == RangeOperation::r_sum){
                return rrr_partial_sums_vector<size_type>(ms_path).noindex(from_idx, to_idx, algo);
            } else {
                return rrr_partial_max_vector<size_type>(ms_path).noindex(from_idx, to_idx, algo);
            }
        }
        static void rrr_nonidex_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage, const RangeAlgorithm algo, const RangeOperation op){

            if(op == RangeOperation::r_sum){
                _noindex_profile<rrr_partial_sums_vector<size_type>>(ms_path, nqueries, range_size, from_idx_max, time_usage, algo);
            } else {
                _noindex_profile<rrr_partial_max_vector<size_type>>(ms_path, nqueries, range_size, from_idx_max, time_usage, algo);
            }
        }

        template <typename op_cls, typename idx_cls>
        static void _indexed_profile(op_cls& pmax, idx_cls& rmq, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                const int block_size, counter_t& time_usage, const RangeAlgorithm algo){

            auto comp_start = timer::now();
            for (int k = 0; k < nqueries; k++) {
                size_type start = random_index(from_idx_max);
                pmax.indexed(rmq, start, start + range_size, (size_type) block_size, algo);
            }
            time_usage.register_now("algorithm", comp_start);
            add_times_to(time_usage, pmax.m_time_usage);
        }

        static rqres_t rrr_indexed(const string ms_path, const string ridx_path,
                const size_type from_idx, const size_type to_idx, const int block_size,
                const RangeAlgorithm algo, const RangeOperation op){

            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);

            if(op == RangeOperation::r_max){
                sdsl::rmq_succinct_sct<false> rmq(&ridx);
                return rrr_partial_max_vector<size_type>(ms_path).indexed(rmq, from_idx, to_idx, (size_type) block_size, algo);
            } else {
                return rrr_partial_sums_vector<size_type>(ms_path).indexed(ridx, from_idx, to_idx, (size_type) block_size, algo);
            }
        }
        static void rrr_indexed_profile(const string ms_path, const string ridx_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                const int block_size, counter_t& time_usage, const RangeAlgorithm algo, const RangeOperation op){

            sdsl::int_vector<64> ridx;
            _load_time_ridx(ridx, ridx_path, time_usage);

            if(op == RangeOperation::r_max){
                rrr_partial_max_vector<size_type> pmax(ms_path);
                auto comp_start = timer::now();
                sdsl::rmq_succinct_sct<false> rmq(&ridx);
                time_usage.register_now("rmq_init", comp_start);

                _indexed_profile<rrr_partial_max_vector<size_type>>(pmax, rmq, nqueries, range_size, from_idx_max, block_size, time_usage, algo);

            } else if (op == RangeOperation::r_sum){
                rrr_partial_sums_vector<size_type> psum(ms_path);
                _indexed_profile<rrr_partial_sums_vector<size_type>>(psum, ridx, nqueries, range_size, from_idx_max, block_size, time_usage, algo);
            } else
                throw string{"Operation max not implemented with index."};
        }

        static rqres_t none_indexed(const string ms_path, const string ridx_path,
                const size_type from_idx, const size_type to_idx, const int block_size,
                const RangeAlgorithm algo, const RangeOperation op){

            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);

            counter_t tusage;
            if(op == RangeOperation::r_max){
                none_partial_max_vector<size_type> pmax(ms_path);
                sdsl::rmq_succinct_sct<false> rmq(&ridx);
                rqres_t answer = pmax.indexed(rmq, from_idx, to_idx, (size_type) block_size, algo);
                return answer;
            } else if (op == RangeOperation::r_sum) {
                none_partial_sums_vector<size_type> psum(ms_path);
                rqres_t answer = psum.indexed(ridx, from_idx, to_idx, (size_type) block_size, algo);
                return answer;
            } else
                throw string{"Operation max not implemented with index."};
        }
        static void none_indexed_profile(const string ms_path, const string ridx_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                const int block_size, counter_t& time_usage, const RangeAlgorithm algo, const RangeOperation op){

            sdsl::int_vector<64> ridx;
            _load_time_ridx(ridx, ridx_path, time_usage);


            if(op == RangeOperation::r_max){
                none_partial_max_vector<size_type> pmax(ms_path);

                auto comp_start = timer::now();
                sdsl::rmq_succinct_sct<false> rmq(&ridx);
                time_usage.register_now("rmq_init", comp_start);
                _indexed_profile<none_partial_max_vector<size_type>>(pmax, rmq, nqueries, range_size, from_idx_max, block_size, time_usage, algo);
            } else if (op == RangeOperation::r_sum) {
                none_partial_sums_vector<size_type> psum(ms_path);
                _indexed_profile<none_partial_sums_vector<size_type>>(psum, ridx, nqueries, range_size, from_idx_max, block_size, time_usage, algo);
            } else
                throw string{"Operation max not implemented with index."};
        }
    };


    template<typename vec_type, typename it_type>
    class rle_rq_dispatcher : rq_dispatcher {
        template<typename opname>
        static void __no_ridx_op_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                const RangeAlgorithm algo,
                counter_t& time_usage){

            auto comp_start = timer::now();
            std::ifstream in{ms_path, std::ios::binary};
            vec_type ms(in);
            it_type* it = new it_type(ms);
            time_usage.register_now("load_ms", comp_start);

            time_usage.register_now("load_partial_ops", timer::now());

            comp_start = timer::now();
            opname p_op(ms, it);
            time_usage.register_now("pop_init", comp_start);

            comp_start = timer::now();
            for (int k = 0; k < nqueries; k++){
                size_type start_idx = random_index(from_idx_max);
                size_type end_idx = start_idx + range_size;
                p_op.noindex(start_idx, end_idx, algo);
            }
            time_usage.register_now("algorithm", comp_start);
        }

        template<typename opname>
        static void __ridx_op_profile(const string ms_path, const string ridx_path,
                const size_type nqueries, const size_type range_size, const size_type from_idx_max,
                const size_type block_size,
                counter_t& time_usage){

            auto comp_start = timer::now();
            std::ifstream in{ms_path, std::ios::binary};
            vec_type ms(in);
            it_type* it = new it_type(ms);
            time_usage.register_now("load_ms", comp_start);

            comp_start = timer::now();
            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);
            time_usage.register_now("load_partial_ops", comp_start);

            comp_start = timer::now();
            opname p_op(ms, it);
            time_usage.register_now("pop_init", comp_start);

            comp_start = timer::now();
            for (int k = 0; k < nqueries; k++) {
                size_type start_idx = random_index(from_idx_max);
                size_type end_idx = start_idx + range_size;
                p_op.indexed(ridx, start_idx, end_idx, block_size);
            }
            time_usage.register_now("algorithm", comp_start);
        }

    public:
        static rqres_t noindex(const string ms_path,
                const size_type from_idx, const size_type to_idx,
                const RangeAlgorithm algo, const RangeOperation op){
            std::ifstream in{ms_path, std::ios::binary};
            vec_type ms(in);
            it_type* it = new it_type(ms);

            if(op == RangeOperation::r_sum){
                rle_partial_sums_vector<vec_type, it_type, size_type> psum(ms, it);
                return psum.noindex(from_idx, to_idx, algo);
            } else {
                rle_partial_max_vector<vec_type, it_type, size_type> pmax(ms, it);
                return pmax.noindex(from_idx, to_idx, algo);
            }
        }
        static void noindex_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage, const RangeAlgorithm algo, const RangeOperation op){

            if(op == RangeOperation::r_sum){
                __no_ridx_op_profile<rle_partial_sums_vector<vec_type, it_type, size_type>>(ms_path, nqueries, range_size, from_idx_max, algo, time_usage);
            } else {
                __no_ridx_op_profile<rle_partial_max_vector<vec_type, it_type, size_type>>(ms_path, nqueries, range_size, from_idx_max, algo, time_usage);
            }
        }

        static rqres_t indexed(const string ms_path, const string ridx_path,
                const size_type from_idx, const size_type to_idx, const size_type block_size,
                const RangeAlgorithm algo, const RangeOperation op){
            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);
            std::ifstream in{ms_path, std::ios::binary};
            vec_type ms(in);
            it_type* it = new it_type(ms);

            if(op == RangeOperation::r_sum){
                rle_partial_sums_vector<vec_type, it_type, size_type> psum(ms, it);
                rqres_t answer = psum.indexed(ridx, from_idx, to_idx, (size_type) block_size);
                return answer;
            } else {
                rle_partial_max_vector<vec_type, it_type, size_type> pmax(ms, it);
                rqres_t answer = pmax.indexed(ridx, from_idx, to_idx, (size_type) block_size);
                return answer;
            }
        }
        static void indexed_profile(const string ms_path, const string ridx_path,
                const size_type nqueries, const size_type range_size, const size_type from_idx_max,
                const size_type block_size,
                counter_t& time_usage, const RangeAlgorithm algo, const RangeOperation op){

            if(op == RangeOperation::r_sum){
                __ridx_op_profile<rle_partial_sums_vector<vec_type, it_type, size_type>>(ms_path, ridx_path, nqueries, range_size, from_idx_max, block_size, time_usage);
            } else {
                __ridx_op_profile<rle_partial_max_vector<vec_type, it_type, size_type>>(ms_path, ridx_path, nqueries, range_size, from_idx_max, block_size, time_usage);
            }
        }
    };
}
