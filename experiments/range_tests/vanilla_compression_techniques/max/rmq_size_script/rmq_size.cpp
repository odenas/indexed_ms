#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <future>
#include <thread>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/partial_sums_vector.hpp"
#include "fd_ms/partial_max_vector.hpp"
#include "fd_ms/help.hpp"
#include "../malloc_count/malloc_count.h"

using namespace std;
using namespace fdms;

typedef typename StreeOhleb<>::size_type size_type;
typedef sdsl::int_vector_buffer<1> buff_vec_t;


class InputFlags {
public:
    bool time_usage;
    size_t block_size;

    InputFlags() { }

    InputFlags(const InputFlags& f) : time_usage{f.time_usage}, block_size{f.block_size} { }

    InputFlags(OptParser input) :
    time_usage{input.getCmdOption("-time_usage") == "1"}
    {
        string bs{input.getCmdOption("-block_size")};
        try {
            block_size = static_cast<size_t> (std::stoi(bs));
        } catch (string m) {
            cerr << ("Error: " + m +
                    "A positive block size is needed (-block_size ). ");
            exit(1);
        }
    }
};


size_type sum(const sdsl::int_vector<64>& ridx){
    size_type a = 0;
    for(size_type i = 0; i < ridx.size() - 2; i += 2)
        a += ridx[i];
    cerr << ridx.size() << " -> " << a << endl;
    return a;
}

void comp1(const string& ridx_path, const size_t block_size){
    sdsl::int_vector<64> ridx;
    sdsl::memory_monitor::start();

    sdsl::load_from_file(ridx, ridx_path);
    cerr << ridx.size() << endl;
    (cout << "rdx," << sdsl::memory_monitor::peak() << "," << block_size << endl);

    sdsl::rmq_succinct_sct<false> rmq(&ridx);
    cerr << rmq(1, 2) << endl;
    (cout << "rmq," << sdsl::memory_monitor::peak() << "," << block_size << endl);
    sdsl::memory_monitor::stop();
}


size_type abs_point() {
    malloc_count_reset_peak();
    return (size_type) malloc_count_peak();
}

size_type diff_from(const size_type from){
    size_type to = abs_point();
    if (from > to)
        throw string{"from peak (" + to_string(from) + ") < to peak (" + to_string(to) + ")"};
    return (size_type) (to - from);
}

void comp2(const string& ridx_path, const size_t block_size){
    sdsl::int_vector<64> ridx;

    size_type from = abs_point();
    sdsl::load_from_file(ridx, ridx_path);
    cerr << ridx.size() << endl;
    (cout << "rdx," << diff_from(from) << "," << block_size << "," << endl);

    from = abs_point();
    sdsl::rmq_succinct_sct<false> rmq(&ridx);
    cerr << rmq(1, 2) << endl;
    (cout << "rmq," << diff_from(from) << "," << block_size << endl);
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    string ms_path;
    string ridx_path;
    InputFlags flags;

    if (argc == 1) {
        (cerr << "Compute block partial sums of bit vector vector (see also 'range_queries.x').\n"
                << "Creates file <ms_path>.<block_size>.ridx in the dir of <ms_path>\n"
                << "Args:\n"
                //<< help__ms_path
                << "\t-ridx_path\n"
                << "\t-block_size <positive int>: the block size; smaller values ~ larger index & faster range query results.\n"
                << endl);
        exit(0);
    } else {
        ridx_path = input.getCmdOption("-ridx_path");
        flags = InputFlags(input);
    }
    cout << "itm,size,bsize" << endl;
    //comp1(ridx_path, flags.block_size);
    comp2(ridx_path, flags.block_size);

}
