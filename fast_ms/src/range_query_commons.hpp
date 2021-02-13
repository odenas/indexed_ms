#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/select_support.hpp>

#include "fd_ms/counter.hpp"
#include "fd_ms/partial_sums_vector.hpp"

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
    };

    class sdsl_rq_dispatcher : rq_dispatcher {
        static void _load_time_ridx(sdsl::int_vector<64>& ridx, const string ridx_path, counter_t& time_usage) {
            auto ds_start = timer::now();
            sdsl::load_from_file(ridx, ridx_path);
            time_usage.register_now("load_partial_sums", ds_start);
        }

    public:
        static size_type none_noindex(const string ms_path,
                const size_type from_idx, const size_type to_idx,
                const RangeAlgorithm algo, const RangeOperation op){
            counter_t  time_usage;
            size_type answer = 0, check_answer = 0;
            if(op == RangeOperation::r_sum){
                return none_partial_sums_vector<size_type>(ms_path, time_usage).noindex(from_idx, to_idx, algo);
            } else {
                return none_partial_max_vector<size_type>(ms_path).noindex(from_idx, to_idx, algo);
            }
        }
        static void none_noindex_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage, const RangeAlgorithm algo, const RangeOperation op){
            if(op == RangeOperation::r_sum){
                none_partial_sums_vector<size_type> psum(ms_path, time_usage);

                auto comp_start = timer::now();
                for (int k = 0; k < nqueries; k++) {
                    size_type start = random_index(from_idx_max);
                    psum.noindex(start, start + range_size, algo);
                }
                time_usage.register_now("algorithm", comp_start);
            } else {
                none_partial_max_vector<size_type> pmax(ms_path, time_usage);

                auto comp_start = timer::now();
                for (int k = 0; k < nqueries; k++) {
                    size_type start = random_index(from_idx_max);
                    pmax.noindex(start, start + range_size, algo);
                }
                time_usage.register_now("algorithm", comp_start);
            }
        }

        static size_type rrr_noindex(const string ms_path,
                const size_type from_idx, const size_type to_idx,
                const RangeAlgorithm algo, const RangeOperation op){
            counter_t  time_usage;
            if(op == RangeOperation::r_sum){
                rrr_partial_sums_vector<size_type> psum(ms_path, time_usage);
                size_type answer = psum.noindex_range_sum(from_idx, to_idx, algo);
                return answer;
            } else {
                rrr_partial_max_vector<size_type> pmax(ms_path);
                size_type answer = pmax.noindex(from_idx, to_idx, algo);
                return answer;
            }
        }
        static void rrr_nonidex_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage, const RangeAlgorithm algo, const RangeOperation op){

            if(op == RangeOperation::r_sum){
                rrr_partial_sums_vector<size_type> psum(ms_path, time_usage);

                auto comp_start = timer::now();
                for (int k = 0; k < nqueries; k++) {
                    size_type start = random_index(from_idx_max);
                    psum.noindex_range_sum(start, start + range_size, algo);
                }
                time_usage.register_now("algorithm", comp_start);
            } else {
                rrr_partial_max_vector<size_type> pmax(ms_path, time_usage);

                auto comp_start = timer::now();
                for (int k = 0; k < nqueries; k++) {
                    size_type start = random_index(from_idx_max);
                    pmax.noindex(start, start + range_size, algo);
                }
                time_usage.register_now("algorithm", comp_start);
            }
        }

        static size_type rrr_indexed(const string ms_path, const string ridx_path,
                const size_type from_idx, const size_type to_idx, const int block_size,
                const RangeAlgorithm algo, const RangeOperation op){

            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);
            if(op == RangeOperation::r_max){
                rrr_partial_max_vector<size_type> pmax(ms_path);
                sdsl::rmq_succinct_sct<false> rmq(&ridx);
                sdsl::rrr_vector<>::rank_1_type rb(&pmax.m_ms);
                counter_t tusage;
                size_type answer = pmax.indexed(rmq, rb, from_idx, to_idx, (size_type) block_size, algo, tusage);
                return answer;
            } else {
                return rrr_partial_sums_vector<size_type>(ms_path).indexed_range_sum(ridx, from_idx, to_idx, (size_type) block_size, algo);
            }
        }

        static void rrr_indexed_profile(const string ms_path, const string ridx_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                const int block_size, counter_t& time_usage, const RangeAlgorithm algo, const RangeOperation op){

            sdsl::int_vector<64> ridx;
            _load_time_ridx(ridx, ridx_path, time_usage);

            if(op == RangeOperation::r_max){
                auto comp_start = timer::now();
                rrr_partial_max_vector<size_type> pmax(ms_path, time_usage);
                sdsl::rmq_succinct_sct<false> rmq(&ridx);
                sdsl::rrr_vector<>::rank_1_type rb(&pmax.m_ms);
                time_usage.register_now("rmq_and_rank_init", comp_start);

                time_usage.reg["bit_range"] = static_cast<size_type>(0);
                comp_start = timer::now();
                for (int k = 0; k < nqueries; k++) {
                    size_type start = random_index(from_idx_max);
                    pmax.indexed(rmq, rb, start, start + range_size, (size_type) block_size, algo, time_usage);
                }
                time_usage.register_now("algorithm", comp_start);
            } else if (op == RangeOperation::r_sum){
                rrr_partial_sums_vector<size_type> psum(ms_path, time_usage);
                auto comp_start = timer::now();
                for (int k = 0; k < nqueries; k++) {
                    size_type start = random_index(from_idx_max);
                    psum.indexed_range_sum(ridx, start, start + range_size, (size_type) block_size, algo);
                }
                time_usage.register_now("algorithm", comp_start);
            } else
                throw string{"Operation max not implemented with index."};
        }

        static void none_indexed_profile(const string ms_path, const string ridx_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                const int block_size, counter_t& time_usage, const RangeAlgorithm algo, const RangeOperation op){

            sdsl::int_vector<64> ridx;
            _load_time_ridx(ridx, ridx_path, time_usage);

            if(op == RangeOperation::r_max){
                auto comp_start = timer::now();
                none_partial_max_vector<size_type> pmax(ms_path);
                sdsl::rmq_succinct_sct<false> rmq(&ridx);
                sdsl::bit_vector::rank_1_type rb(&pmax.m_ms);
                time_usage.register_now("rmq_and_rank_init", comp_start);

                std::vector<string> _keys = {
                    "algorithm.trivial_case",
                    "algorithm.rmq_scan", "algorithm.rmq_query",
                    "algorithm.trivial_scan",
                    "algorithm.trivial_scan.1", "algorithm.trivial_scan.2"
                };
                for(auto k: _keys)
                    time_usage.register_now(k, timer::now());

                for(auto k: {"range.int", "range.bit", "range.block", "range.i_block"})
                    time_usage.reg[k] = static_cast<size_type>(0);

                comp_start = timer::now();
                for (int k = 0; k < nqueries; k++) {
                    size_type start = random_index(from_idx_max);
                    pmax.indexed(rmq, rb, start, start + range_size, (size_type) block_size, algo, time_usage);
                }
                time_usage.register_now("algorithm", comp_start);
            } else if (op == RangeOperation::r_sum) {
                none_partial_sums_vector<size_type> psum(ms_path, time_usage);
                auto comp_start = timer::now();
                for (int k = 0; k < nqueries; k++) {
                    size_type start = random_index(from_idx_max);
                    psum.indexed(ridx, start, start + range_size, (size_type) block_size, algo);
                }
                time_usage.register_now("algorithm", comp_start);
            } else
                throw string{"Operation max not implemented with index."};
        }
        static size_type none_indexed(const string ms_path, const string ridx_path,
                const size_type from_idx, const size_type to_idx, const int block_size,
                const RangeAlgorithm algo, const RangeOperation op){

            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);

            if(op == RangeOperation::r_max){
                none_partial_max_vector<size_type> pmax(ms_path);
                sdsl::rmq_succinct_sct<false> rmq(&ridx);
                sdsl::bit_vector::rank_1_type rb(&pmax.m_ms);
                counter_t tusage;
                size_type answer = pmax.indexed(rmq, rb, from_idx, to_idx, (size_type) block_size, algo, tusage);
                return answer;
            } else {
                none_partial_sums_vector<size_type> psum(ms_path);
                size_type answer = psum.indexed(ridx, from_idx, to_idx, (size_type) block_size, algo);
                return answer;
            }
        }
    };


    template<typename vec_type, typename it_type>
    class rle_rq_dispatcher : rq_dispatcher {
        template<typename opname>
        static void __no_ridx_op_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                const RangeAlgorithm algo,
                counter_t& time_usage){

            time_usage.register_now("init.load_psums", timer::now());

            auto comp_start = timer::now();
            std::ifstream in{ms_path, std::ios::binary};
            vec_type ms(in);
            it_type* it = new it_type(ms);
            time_usage.register_now("load_ms", comp_start);

            time_usage.register_now("load_partial_sums", timer::now());

            comp_start = timer::now();
            opname p_op(ms, it);
            for (int k = 0; k < nqueries; k++){
                size_type start_idx = random_index(from_idx_max);
                size_type end_idx = start_idx + range_size;
                p_op.noindex(start_idx, end_idx, algo);
            }
            time_usage.register_now("algorithm", comp_start);
        }

    public:
        static size_type noindex(const string ms_path,
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

        static void indexed_profile(const string ms_path, const string ridx_path,
                const size_type nqueries, const size_type range_size, const size_type from_idx_max,
                const size_type block_size,
                counter_t& time_usage){

            auto comp_start = timer::now();
            std::ifstream in{ms_path, std::ios::binary};
            vec_type ms(in);
            time_usage.register_now("load_ms", comp_start);

            auto ds_start = timer::now();
            it_type* it = new it_type(ms);
            time_usage.register_now("select_init", ds_start);

            ds_start = timer::now();
            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);
            time_usage.register_now("init.load_psums", ds_start);

            ds_start = timer::now();
            rle_partial_sums_vector<vec_type, it_type, size_type> psum(ms, it);
            time_usage.register_now("pmax_init", ds_start);

            comp_start = timer::now();
            for (int k = 0; k < nqueries; k++) {
                size_type start_idx = random_index(from_idx_max);
                size_type end_idx = start_idx + range_size;
                psum.indexed(ridx, start_idx, end_idx, block_size);
            }
            time_usage.register_now("algorithm", comp_start);
        }

        static size_type indexed(const string ms_path, const string ridx_path,
                const size_type from_idx, const size_type to_idx, const size_type block_size,
                const RangeAlgorithm algo, const RangeOperation op){
            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);
            std::ifstream in{ms_path, std::ios::binary};
            vec_type ms(in);
            it_type* it = new it_type(ms);

            if(op == RangeOperation::r_sum){
                rle_partial_sums_vector<vec_type, it_type, size_type> psum(ms, it);
                size_type answer = psum.indexed(ridx, from_idx, to_idx, (size_type) block_size);
                return answer;
            } else {
                rle_partial_max_vector<vec_type, it_type, size_type> pmax(ms, it);
                size_type answer = pmax.indexed(ridx, from_idx, to_idx, (size_type) block_size);
                return answer;
            }
        }
    };
}
