#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <future>
#include <thread>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/counter.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/partial_sums_vector.hpp"
#include "fd_ms/help.hpp"

using namespace std;
using namespace fdms;

typedef typename StreeOhleb<>::size_type size_type;
typedef Counter<size_type> counter_t;
typedef sdsl::int_vector_buffer<1> buff_vec_t;


class InputFlags {
private:

    void check() const {
        if (block_size < 0) {
            cerr << "positive block size needed. got: " << block_size << endl;
            exit(1);
        }
    }

public:
    bool time_usage; // can be max/min etc.
    size_t block_size;

    InputFlags() { }

    InputFlags(const InputFlags& f) : time_usage{f.time_usage}, block_size{f.block_size} { }

    InputFlags(bool time_, size_t block_size) : time_usage{time_}, block_size{block_size}
    {
        check();
    }

    InputFlags(OptParser input) : time_usage{input.getCmdOption("-time_usage") == "1"}
    {
        string bs{input.getCmdOption("-block_size")};
        try {
            block_size = static_cast<size_t> (std::stoi(bs));
        } catch (string m) {
            cerr << ("Error: " + m +
                    "A positive block size is needed (-block_size ). ");
            exit(1);
        }
        check();
    }
};

static void show_MS(const InputSpec ispec, std::ostream& out) {
    buff_vec_t ms(ispec.ms_fname, std::ios::in);
    //sdsl::bit_vector ms;
    //sdsl::load_from_file(ms, ispec.ms_fname);

    size_type k = 0;
    for (size_type i = 0; i < ms.size(); i++) {
        if (ms[i] == 1) {
            out << i - (2 * k) << " ";
            k += 1;
        }
    }
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    string ms_path;
    InputFlags flags;
    counter_t time_usage{};

    if (argc == 1) {
        (cerr << "Dump an index needed to speedup the range queries (see 'range_queries.x').\n"
                << "Creates file <ms_path>.<block_size>.ridx in the dir of <ms_path>\n"
                << "Args:\n"
                << help__ms_path
                << "\t-block_size <positive int>: the block size; smaller values ~ larger index & faster range query results.\n"
                << help__time_usage
                << endl);
        exit(0);

        ms_path = {"/home/brt/code/matching_statistics/indexed_ms/fast_ms/tests/"};
        flags = InputFlags(
                false, // time
                4 // block_size
                );
    } else {
        ms_path = input.getCmdOption("-ms_path");
        flags = InputFlags(input);
    }

    //cout << "ms_path: " << ms_path << endl;

    auto comp_start = timer::now();
    try{
        partial_sums_vector<size_type, sdsl::bit_vector, sdsl::bit_vector::select_1_type>::dump(ms_path, flags.block_size);
    } catch (string s) {
        cerr << "Couldn't dump ms partials sums. Reason: " << endl;
        cerr << s << endl;
    }
    time_usage.register_now("total_time", comp_start);

    if (flags.time_usage) {
        cerr << "dumping reports ..." << endl;
        for (auto item : time_usage.reg)
            (cout << item.first << "," << item.second << endl);
    }
    return 0;
}
