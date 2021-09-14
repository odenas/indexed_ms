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


void comp2(const string& ridx_path, const size_t block_size){
    sdsl::int_vector<64> ridx;
    sdsl::load_from_file(ridx, ridx_path);
    cout << "ridx," << size_in_bytes(ridx) << endl;

    sdsl::rmq_succinct_sct<false> rmq(&ridx);
    rmq(1, 2);
    cout << "rmq," << size_in_bytes(rmq) << endl;
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
                << "\t-ridx_path\n"
                << "\t-block_size <positive int>: the block size; smaller values ~ larger index & faster range query results.\n"
                << endl);
        exit(0);
    } else {
        ridx_path = input.getCmdOption("-ridx_path");
        flags = InputFlags(input);
    }
    cout << "dstructure,size" << endl;
    comp2(ridx_path, flags.block_size);

}
