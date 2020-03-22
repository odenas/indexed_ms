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
        static void _load_time_ridx(sdsl::int_vector<64>& ridx, const string ridx_path, counter_t& time_usage) {
            auto ds_start = timer::now();
            sdsl::load_from_file(ridx, ridx_path);
            time_usage.register_now("load_partial_sums", ds_start);
        }

    public:
        static size_type none_noindex(const string ms_path,
                const size_type from_idx, const size_type to_idx,
                const bool check, const IndexedAlgorithm algo, const RangeOperation op){
            counter_t  time_usage;
            size_type answer = 0, check_answer = 0;
            if(op == RangeOperation::r_sum){
                none_partial_sums_vector<size_type> psum(ms_path, time_usage);

                answer = psum.noindex_range_sum(from_idx, to_idx, algo);
                return (check ? __check_outcome(answer, psum.check_range_sum(from_idx, to_idx)) : answer);
            } else {
                none_partial_max_vector<size_type> pmax(ms_path);

                answer = pmax.noindex_range_max(from_idx, to_idx, algo);
                return (check ? __check_outcome(answer, pmax.check_range_max(from_idx, to_idx)) : answer);
            }
        }
        static void none_noindex_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage, const IndexedAlgorithm algo, const RangeOperation op){
            if(op == RangeOperation::r_sum){
                none_partial_sums_vector<size_type> psum(ms_path, time_usage);

                auto comp_start = timer::now();
                for (int k = 0; k < nqueries; k++) {
                    size_type start = random_index(from_idx_max);
                    psum.noindex_range_sum(start, start + range_size, algo);
                }
                time_usage.register_now("algorithm", comp_start);
            } else {
                none_partial_max_vector<size_type> pmax(ms_path, time_usage);

                auto comp_start = timer::now();
                for (int k = 0; k < nqueries; k++) {
                    size_type start = random_index(from_idx_max);
                    pmax.noindex_range_max(start, start + range_size, algo);
                }
                time_usage.register_now("algorithm", comp_start);
            }
        }

        static size_type rrr_noindex(const string ms_path,
                const size_type from_idx, const size_type to_idx,
                const bool check, const IndexedAlgorithm algo, const RangeOperation op){
            counter_t  time_usage;
            if(op == RangeOperation::r_sum){
                rrr_partial_sums_vector<size_type> psum(ms_path, time_usage);
                size_type answer = psum.noindex_range_sum(from_idx, to_idx, algo);
                return (check ? __check_outcome(answer, psum.check_range_sum(from_idx, to_idx)) : answer);
            } else {
                rrr_partial_max_vector<size_type> pmax(ms_path);
                size_type answer = pmax.noindex_range_max(from_idx, to_idx, algo);
                return (check ? __check_outcome(answer, pmax.check_range_max(from_idx, to_idx)) : answer);
            }
        }
        static void rrr_nonidex_profile(const string ms_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                counter_t& time_usage, const IndexedAlgorithm algo, const RangeOperation op){

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
                    pmax.noindex_range_max(start, start + range_size, algo);
                }
                time_usage.register_now("algorithm", comp_start);
            }
        }

        static size_type rrr_indexed(const string ms_path, const string ridx_path,
                const size_type from_idx, const size_type to_idx, const int block_size,
                const bool check, const IndexedAlgorithm algo, const RangeOperation op){
            counter_t  time_usage;
            if(op == RangeOperation::r_max)
                throw string{"Operation max not implemented with index."};

            rrr_partial_sums_vector<size_type> psum(ms_path, time_usage);

            sdsl::int_vector<64> ridx;
            sdsl::load_from_file(ridx, ridx_path);
            size_type answer = psum.indexed_range_sum(ridx, from_idx, to_idx, (size_type) block_size, algo);
            return (check ? __check_outcome(answer, psum.check_range_sum(from_idx, to_idx)) : answer);
        }
        static void rrr_indexed_profile(const string ms_path, const string ridx_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                const int block_size, counter_t& time_usage, const IndexedAlgorithm algo, const RangeOperation op){
            if(op == RangeOperation::r_max)
                throw string{"Operation max not implemented with index."};

            rrr_partial_sums_vector<size_type> psum(ms_path, time_usage);

            sdsl::int_vector<64> ridx;
            _load_time_ridx(ridx, ridx_path, time_usage);

            auto comp_start = timer::now();
            for (int k = 0; k < nqueries; k++) {
                size_type start = random_index(from_idx_max);
                psum.indexed_range_sum(ridx, start, start + range_size, (size_type) block_size, algo);
            }
            time_usage.register_now("algorithm", comp_start);
        }

        static void none_indexed_profile(const string ms_path, const string ridx_path, const size_type nqueries,
                const size_type range_size, const size_type from_idx_max,
                const int block_size, counter_t& time_usage, const IndexedAlgorithm algo, const RangeOperation op){
            if(op == RangeOperation::r_max)
                throw string{"Operation max not implemented with index."};


            none_partial_sums_vector<size_type> psum(ms_path, time_usage);

            sdsl::int_vector<64> ridx;
            _load_time_ridx(ridx, ridx_path, time_usage);

            auto comp_start = timer::now();
            for (int k = 0; k < nqueries; k++) {
                size_type start = random_index(from_idx_max);
                psum.indexed_range_sum(ridx, start, start + range_size, (size_type) block_size, algo);
            }
            time_usage.register_now("algorithm", comp_start);
        }
        static size_type none_indexed(const string ms_path, const string ridx_path,
                const size_type from_idx, const size_type to_idx, const int block_size,
                const bool check, const IndexedAlgorithm algo, const RangeOperation op){
            counter_t  time_usage;
            if(op == RangeOperation::r_max){
                none_partial_max_vector<size_type> pmax(ms_path, time_usage);

                sdsl::int_vector<64> ridx;
                sdsl::load_from_file(ridx, ridx_path);
                size_type answer = pmax.indexed_range_max(ridx, from_idx, to_idx, (size_type) block_size, algo);
                return (check ? __check_outcome(answer, pmax.check_range_max(from_idx, to_idx)) : answer);
            } else {
                none_partial_sums_vector<size_type> psum(ms_path, time_usage);

                sdsl::int_vector<64> ridx;
                sdsl::load_from_file(ridx, ridx_path);
                size_type answer = psum.indexed_range_sum(ridx, from_idx, to_idx, (size_type) block_size, algo);
                return (check ? __check_outcome(answer, psum.check_range_sum(from_idx, to_idx)) : answer);
            }
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