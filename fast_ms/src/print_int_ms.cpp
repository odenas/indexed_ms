/*
print sections of an MS bit vector
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "sdsl/vectors.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/util.hpp"

#include "fd_ms/help.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/counter.hpp"

using namespace fdms;
using namespace std;

typedef unsigned long long size_type;
typedef sdsl::int_vector_buffer<1> buff_vec_t;

int comp(const string ms_path) {
    buff_vec_t ms(ms_path, std::ios::in);
    //sdsl::bit_vector ms;
    //sdsl::load_from_file(ms, ms_path);

    size_type k = 0, max_k = ms.size() / 2;
    for (size_type i = 0; i < ms.size(); i++) {
        if (ms[i] == 1) {
            cout << i - (2 * k);
            k += 1;
            if(k < max_k)
                cout << " ";
        }
    }
    cout << endl;
    return 0;
}

int main(int argc, char **argv) {
    if(argc == 1){
        (cerr << "Print the given ms bit-vector\n"
              << "Args:\n"
              << help__ms_path
              << endl);
        exit(0);
    }

    return comp(OptParser(argc, argv).getCmdOption("-ms_path"));
}
