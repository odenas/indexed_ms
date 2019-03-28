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
#include "fd_ms/counter.hpp"

#include "rlcsa/bits/bitvector.h"
#include "rlcsa/bits/rlevector.h"
#include "rlcsa/bits/deltavector.h"
#include "rlcsa/bits/succinctvector.h"
#include "rlcsa/bits/nibblevector.h"


using namespace fdms;
using namespace std;
using timer = std::chrono::high_resolution_clock;
typedef unsigned long long size_type;
typedef std::map<size_type, size_type> histo_t;

template<typename enc_type>
void comp(const string ms_path, const string out_suffix){
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);
    //sdsl::rrr_vector<> c_ms(ms);
    //sdsl::hyb_vector<> c_ms(ms);
    enc_type c_ms(ms);
    sdsl::store_to_file(c_ms, ms_path + out_suffix);
}

template<typename enc_type>
size_type fill_encoder(sdsl::bit_vector ms, enc_type& encoder, histo_t& freq){
    size_type nz = 0, i = 0, n_runs = 0;
    while(i < ms.size()){
        if(ms[i] == 0)
            nz += 1;
        else{
            assert (i >= nz);
            if (nz > 0){
                n_runs += 1;
                encoder.addRun(i - nz, nz);
                freq[nz] += 1;
            }
            nz = 0;
        }
        i += 1;
    }
    encoder.flush();
    return n_runs;
}

template<typename vec_type, typename enc_type>
void comp1(const string ms_path, const string suffix){
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);

    enc_type encoder(32);
    histo_t counter;
    size_type n_runs = fill_encoder<enc_type>(ms, encoder, counter);
    vec_type c_ms(encoder, ms.size());
    std::ofstream out{ms_path + "." + suffix, std::ios::binary};
    c_ms.writeTo(out);

    (cerr << n_runs << " runs over "
          << c_ms.getSize() << " elements ("
          << c_ms.getSize() / static_cast<float>(n_runs) << " elements / run)"
          << endl);
    for(auto item : counter)
        cout << item.first << "," << item.second << endl;
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
    comp<sdsl::hyb_vector<>>(ms_path, ".hyb");
    comp<sdsl::rrr_vector<>>(ms_path, ".rrr");
    //comp1<CSA::RLEVector, CSA::RLEEncoder>(ms_path, "rle");
    //comp1<CSA::DeltaVector, CSA::DeltaEncoder>(ms_path, "delta");
    //comp1<CSA::SuccinctVector, CSA::SuccinctEncoder>(ms_path, "succinct");
    //comp1<CSA::NibbleVector, CSA::NibbleEncoder>(ms_path, "nibble");
    return 0;
}

