/*
 compress a freq vector
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "sdsl/bit_vectors.hpp"
#include "sdsl/util.hpp"

#include "fd_ms/help.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/stree_sct3.hpp"

#include "rlcsa/bits/bitvector.h"
#include "rlcsa/bits/rlevector.h"

#include "../malloc_count/malloc_count.h"

using namespace fdms;
using namespace std;

typedef StreeOhleb<> cst_t;
typedef typename cst_t::size_type size_type;
typedef typename CSA::RLEEncoder enc_type;
typedef typename CSA::RLEVector vec_type;


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

size_type fill_encoder(sdsl::int_vector<64> ms, enc_type& encoder){
    size_type no = 0, i = 0, n_runs = 0;
    while(i < ms.size()){
        if(ms[i] == 1)
            no += 1;
        else{
            assert (i >= no);
            if (no > 0){
                n_runs += 1;
                encoder.addRun(i - no, no);
                no = 0;
            }
        }
        i += 1;
    }
    if(no > 0){
        encoder.addRun(i - no, no);
        n_runs += 1;
    }
    encoder.flush();
    return n_runs;
}

size_type comp1(const string ms_path){
    sdsl::int_vector<64> freq;
    sdsl::load_from_file(freq, ms_path);

    size_type from = abs_point();
    CSA::RLEEncoder encoder(32);
    size_type n_runs = fill_encoder(freq, encoder);
    CSA::RLEVector c_ms(encoder, freq.size());
    std::ofstream out{ms_path +  ".rle", std::ios::binary};
    c_ms.writeTo(out);

    (cerr << n_runs << " runs over "
          << c_ms.getSize() << " elements ("
          << c_ms.getSize() / static_cast<float>(n_runs) << " elements / run)"
          << endl);
    return diff_from(from);
}


int main(int argc, char **argv){
    OptParser input(argc, argv);

    if(argc == 1){
        (cerr << "Compress a freq vector with rle-encoding. Creates files <path>.rle\n"
              << "Args:\n"
              << help__ms_path
              << endl);
        exit(0);
    }
    cout << comp1(input.getCmdOption("-ms_path")) << endl;
    return 0;
}
