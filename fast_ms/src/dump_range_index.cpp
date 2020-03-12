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
private:
    static RangeOperation parse_operation(const string c_str){
        std::map<RangeOperation, string> a2s = {
            {RangeOperation::r_sum, "sum"},
            {RangeOperation::r_max, "max"}
        };

        if(c_str == "0" or c_str == "sum")
            return RangeOperation::r_sum;
        for(auto item: a2s){
            if(item.second == c_str)
                return item.first;
        }
        throw string{"bad operation string: " + c_str};
    }

    void check() const {
        if (block_size < 0) {
            cerr << "positive block size needed. got: " << block_size << endl;
            exit(1);
        }
    }

public:
    bool time_usage; // can be max/min etc.
    size_t block_size;
    RangeOperation op;

    InputFlags() { }

    InputFlags(const InputFlags& f) : time_usage{f.time_usage}, block_size{f.block_size}, op{f.op} { }

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
        op = parse_operation(input.getCmdOption("-op"));
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

    if (argc == 1) {
        (cerr << "Compute block partial sums of bit vector vector (see also 'range_queries.x').\n"
                << "Creates file <ms_path>.<block_size>.ridx in the dir of <ms_path>\n"
                << "Args:\n"
                << help__ms_path
                << "\t-block_size <positive int>: the block size; smaller values ~ larger index & faster range query results.\n"
                << help__rangeop
                << endl);
        exit(0);
    } else {
        ms_path = input.getCmdOption("-ms_path");
        flags = InputFlags(input);
    }

    auto comp_start = timer::now();
    try{
        if(flags.op == RangeOperation::r_sum)
            sdsl_partial_sums_vector<sdsl::bit_vector, sdsl::bit_vector::select_1_type, size_type>::dump(ms_path, flags.block_size);
        else if (flags.op == RangeOperation::r_max)
            sdsl_partial_max_vector<sdsl::bit_vector, sdsl::bit_vector::select_1_type, size_type>::dump(ms_path, flags.block_size);
    } catch (string s) {
        cerr << "Couldn't dump ms partials sums. Reason: " << endl;
        cerr << s << endl;
    }
    return 0;
}
