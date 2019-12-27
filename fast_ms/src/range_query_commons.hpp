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

        static size_type __check_outcome(const size_type answer, size_type answer_check){
            if (answer != answer_check)
                throw string{"answer " + to_string(answer) + " != expected answer " + to_string(answer_check)};
            return answer;
        }

    public:
        typedef Counter<size_type> counter_t;
    };

    class sdsl_rq_dispatcher : rq_dispatcher {
    public:
        template<typename ms_type>
        static void trivial_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage){

            typedef typename ms_type::select_1_type ms_sel_1_type;
            sdsl_partial_sums_vector<ms_type, ms_sel_1_type, size_type> psum(ms_path, time_usage);

            auto comp_start = timer::now();
            for (int k = 0; k < nqueries; k++) {
                size_type start = random_index(from_idx_max);
                psum.trivial_range_sum(start, start + range_size);
            }
            time_usage.register_now("algorithm", comp_start);
        }

        template<typename ms_type>
        static size_type trivial(const string ms_path,
                const size_type from_idx, const size_type to_idx, const bool check){
            counter_t  time_usage;
            typedef typename ms_type::select_1_type ms_sel_1_type;
            typedef sdsl_partial_sums_vector<ms_type, ms_sel_1_type, size_type> psum_t;


            size_type answer = psum_t(ms_path, time_usage).trivial_range_sum(from_idx, to_idx);
            return (check ? __check_outcome(answer, psum_t(ms_path, time_usage).trivial_range_sum(from_idx, to_idx)) : answer);
        }

        static void rrr_fast_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage){

            rrr_partial_sums_vector<size_type> psum(ms_path, time_usage);

            auto comp_start = timer::now();
            for (int k = 0; k < nqueries; k++) {
                size_type start = random_index(from_idx_max);
                psum.djamal_range_sum(start, start + range_size);
            }
            time_usage.register_now("algorithm", comp_start);
        }

        static size_type rrr_fast(const string ms_path,
                const size_type from_idx, const size_type to_idx, const bool check){
            counter_t  time_usage;
            rrr_partial_sums_vector<size_type> psum(ms_path, time_usage);

            size_type answer = psum.djamal_range_sum(from_idx, to_idx);
            return (check ? __check_outcome(answer, psum.trivial_range_sum(from_idx, to_idx)) : answer);
        }

        static void none_fast_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage){

            none_partial_sums_vector<size_type> psum(ms_path, time_usage);

            auto comp_start = timer::now();
            for (int k = 0; k < nqueries; k++) {
                size_type start = random_index(from_idx_max);
                psum.djamal_range_sum(start, start + range_size);
            }
            time_usage.register_now("algorithm", comp_start);
        }

        static size_type none_fast(const string ms_path,
                const size_type from_idx, const size_type to_idx, const bool check){
            counter_t  time_usage;
            none_partial_sums_vector<size_type> psum(ms_path, time_usage);

            size_type answer = psum.djamal_range_sum(from_idx, to_idx);
            return (check ? __check_outcome(answer, psum.trivial_range_sum(from_idx, to_idx)) : answer);
        }

        template<typename ms_type>
        static void indexed_profile(const string ms_path, const string ridx_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                const int block_size, counter_t& time_usage){

            typedef typename ms_type::select_1_type ms_sel_1_type;
            sdsl_partial_sums_vector<ms_type, ms_sel_1_type, size_type> psum(ms_path, time_usage);

            auto ds_start = timer::now();
            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);
            time_usage.register_now("load_partial_sums", ds_start);

            auto comp_start = timer::now();
            for (int k = 0; k < nqueries; k++) {
                size_type start = random_index(from_idx_max);
                psum.indexed_range_sum(ridx, start, start + range_size, (size_type) block_size);
            }
            time_usage.register_now("algorithm", comp_start);
        }

        template<typename ms_type>
        static size_type indexed(const string ms_path, const string ridx_path,
                const size_type from_idx, const size_type to_idx, const int block_size, const bool check){
            counter_t  time_usage;

            typedef typename ms_type::select_1_type ms_sel_1_type;
            typedef sdsl_partial_sums_vector<ms_type, ms_sel_1_type, size_type> psum_t;

            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);
            size_type answer = psum_t(ms_path, time_usage).indexed_range_sum(ridx, from_idx, to_idx, (size_type) block_size);
            return (check ? __check_outcome(answer, psum_t(ms_path, time_usage).trivial_range_sum(from_idx, to_idx)) : answer);
        }
    };


    template<typename vec_type, typename it_type>
    class rle_rq_dispatcher : rq_dispatcher {
        static void __no_ridx_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                const size_type block_size,
                counter_t& time_usage){

            auto comp_start = timer::now();
            std::ifstream in{ms_path, std::ios::binary};
            vec_type ms(in);
            it_type* it = new it_type(ms);
            time_usage.register_now("load_ms", comp_start);

            time_usage.register_now("load_partial_sums", timer::now());

            comp_start = timer::now();
            rle_partial_sums_vector<vec_type, it_type, size_type> psum(ms, it);
            for (int k = 0; k < nqueries; k++){
                size_type start_idx = random_index(from_idx_max);
                size_type end_idx = start_idx + range_size;
                if (block_size == 0)
                    psum.trivial_range_sum(start_idx, end_idx);
                else if(block_size == -1)
                    psum.rle_range_sum(start_idx, end_idx);
                else
                    throw string{"Bad block size " + to_string(block_size)};
            }
            time_usage.register_now("algorithm", comp_start);
        }

    public:
        static void trivial_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage){
            return __no_ridx_profile(ms_path, nqueries, range_size, from_idx_max, 0, time_usage);
        }

        static size_type trivial(const string ms_path,
                const size_type from_idx, const size_type to_idx, const bool check){
            std::ifstream in{ms_path, std::ios::binary};
            vec_type ms(in);
            it_type* it = new it_type(ms);
            rle_partial_sums_vector<vec_type, it_type, size_type> psum(ms, it);
            size_type answer = psum.trivial_range_sum(from_idx, to_idx);
            return (check ? __check_outcome(answer, psum.trivial_range_sum(from_idx, to_idx)) : answer);
        }

        static void fast_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage){
            return __no_ridx_profile(ms_path, nqueries, range_size, from_idx_max, -1, time_usage);
        }

        static size_type fast(const string ms_path,
                const size_type from_idx, const size_type to_idx, const bool check){
            std::ifstream in{ms_path, std::ios::binary};
            vec_type ms(in);
            it_type* it = new it_type(ms);
            rle_partial_sums_vector<vec_type, it_type, size_type> psum(ms, it);
            size_type answer = psum.rle_range_sum(from_idx, to_idx);
            return (check ? __check_outcome(answer, psum.trivial_range_sum(from_idx, to_idx)) : answer);
        }

        static void indexed_profile(const string ms_path, const string ridx_path,
                const size_type nqueries, const size_type range_size, const size_type from_idx_max, const size_type block_size,
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
            time_usage.register_now("load_partial_sums", ds_start);

            comp_start = timer::now();
            rle_partial_sums_vector<vec_type, it_type, size_type> psum(ms, it);
            for (int k = 0; k < nqueries; k++) {
                size_type start_idx = random_index(from_idx_max);
                size_type end_idx = start_idx + range_size;
                psum.indexed_range_sum(ridx, start_idx, end_idx, block_size);
            }
            time_usage.register_now("algorithm", comp_start);
        }

        static size_type indexed(const string ms_path, const string ridx_path,
                const size_type from_idx, const size_type to_idx, const size_type block_size,
                const bool check){
            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);

            std::ifstream in{ms_path, std::ios::binary};
            vec_type ms(in);
            it_type* it = new it_type(ms);
            rle_partial_sums_vector<vec_type, it_type, size_type> psum(ms, it);

            size_type answer = psum.indexed_range_sum(ridx, from_idx, to_idx, (size_type) block_size);
            return (check ? __check_outcome(answer, psum.trivial_range_sum(from_idx, to_idx)) : answer);
        }

    };
}