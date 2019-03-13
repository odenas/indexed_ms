/*
 compress a ms vector
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "sdsl/bit_vectors.hpp"
#include "sdsl/util.hpp"

#include "fd_ms/help.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"

#include "rlcsa/bits/bitvector.h"
#include "rlcsa/bits/rlevector.h"

using namespace fdms;
using namespace std;
using timer = std::chrono::high_resolution_clock;
typedef unsigned long long size_type;


void comp(const string ms_path){
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);
    sdsl::rrr_vector<> c_ms(ms);
    sdsl::store_to_file(c_ms, ms_path + ".rrr");
}

void comp1(const string ms_path){
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);

    CSA::RLEVector::Encoder encoder(32);
    CSA::pair_type run = make_pair<size_t, size_t>(0, 0);
    size_t nz = 0, i = 0;
    size_type n_runs = 0, min_run_len = ms.size(), max_run_len = 0;
    while(i < ms.size()){
        if(ms[i] == 0)
            nz += 1;
        else{
            assert (i >= nz);
            if (nz > 0){
                min_run_len = min(min_run_len, static_cast<size_type>(nz));
                max_run_len = max(min_run_len, static_cast<size_type>(nz));
                n_runs += 1;
                encoder.addRun(i - nz, nz);
            }
            nz = 0;
        }
        i += 1;
    }
    encoder.flush();
    CSA::RLEVector c_ms(encoder, ms.size());
    (cerr << n_runs << " runs over "
          << c_ms.getSize() << " elements ("
          << c_ms.getSize() / static_cast<float>(n_runs) << " elements / run)"
          << endl);

    std::ofstream out{ms_path + ".rrr", std::ios::binary};
    c_ms.writeTo(out);
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    string ms_path;

    if(argc == 1){
        (cerr << "Compress a ms vector. Creates files <ms_path>.rrr\n"
              << "Args:\n"
              << help__ms_path
              << endl);
        exit(0);
    } else {
        ms_path = input.getCmdOption("-ms_path");
    }
    comp(ms_path);
    return 0;
}

