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
#include "fd_ms/p_ms_vector.hpp"
#include "fd_ms/stree_sct3.hpp"

#include "rlcsa/bits/bitvector.h"
#include "rlcsa/bits/rlevector.h"
#include "rlcsa/bits/deltavector.h"
#include "rlcsa/bits/succinctvector.h"
#include "rlcsa/bits/nibblevector.h"

#include "../malloc_count/malloc_count.h"

using namespace fdms;
using namespace std;

typedef StreeOhleb<> cst_t;
typedef typename cst_t::size_type size_type;
typedef typename ms_compression::compression_types Compression;
typedef std::map<size_type, size_type> histo_t;


class InputFlags {
public:
    Compression compression;

    InputFlags() { }

    InputFlags(const InputFlags& f) : compression{f.compression} { }

    InputFlags(const Compression compression) : compression{compression} { }

    InputFlags(OptParser input) {
        compression = ms_compression::parse_compression(
            input.getCmdOption("-compression")
        );
    }
};


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

template<typename enc_type>
size_type comp(const string ms_path, const InputFlags& flags){
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);

    size_type from = abs_point();
    enc_type c_ms(ms);
    sdsl::store_to_file(c_ms, ms_path + ms_compression::to_str(flags.compression));
    return diff_from(from);
}

template<typename enc_type>
size_type fill_encoder(sdsl::bit_vector ms, enc_type& encoder, histo_t& freq){
    size_type no = 0, i = 0, n_runs = 0;
    while(i < ms.size()){
        if(ms[i] == 1)
            no += 1;
        else{
            assert (i >= no);
            if (no > 0){
                n_runs += 1;
                encoder.addRun(i - no, no);
                freq[no] += 1;
                no = 0;
            }
        }
        i += 1;
    }
    if(no > 0){
        encoder.addRun(i - no, no);
        n_runs += 1;
        freq[no] += 1;
    }
    encoder.flush();
    return n_runs;
}

template<typename vec_type, typename enc_type>
size_type comp1(const string ms_path, const InputFlags& flags){
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);

    size_type from = abs_point();
    enc_type encoder(32);
    histo_t counter;
    size_type n_runs = fill_encoder<enc_type>(ms, encoder, counter);
    vec_type c_ms(encoder, ms.size());
    std::ofstream out{ms_path +  ms_compression::to_str(flags.compression), std::ios::binary};
    c_ms.writeTo(out);

    (cerr << n_runs << " runs over "
          << c_ms.getSize() << " elements ("
          << c_ms.getSize() / static_cast<float>(n_runs) << " elements / run)"
          << endl);
    //for(auto item : counter)
    //    cout << item.first << "," << item.second << endl;
    return diff_from(from);
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    string ms_path;

    if(argc == 1){
        (cerr << "Compress a ms vector. Creates files <ms_path>.xxx\n"
              << "Args:\n"
              << help__ms_path
              << "\t-compression <compression type>: One of: rrr, hybrid, rle, delta, succint, nibble.\n"
              << endl);
        exit(0);
    }
    InputFlags flags;
    try{
        flags = InputFlags(input);
    }
    catch (string s) {
        cerr << s << endl;
        return 1;
    }
    ms_path = input.getCmdOption("-ms_path");
    //size_type mem_mark = abs_point();
    switch(flags.compression)
    {
        case Compression::hybrid:
            cout << comp<sdsl::hyb_vector<>>(ms_path, flags) << endl;
            break;
        case Compression::rrr:
            cout << comp<sdsl::rrr_vector<>>(ms_path, flags) << endl;
            break;
        case Compression::rle:
            cout << comp1<CSA::RLEVector, CSA::RLEEncoder>(ms_path, flags) << endl;
            break;
        case Compression::delta:
            cout << comp1<CSA::DeltaVector, CSA::DeltaEncoder>(ms_path, flags) << endl;
            break;
        case Compression::succint:
            cout << comp1<CSA::SuccinctVector, CSA::SuccinctEncoder>(ms_path, flags) << endl;
            break;
        case Compression::nibble:
            cout << comp1<CSA::NibbleVector, CSA::NibbleEncoder>(ms_path, flags) << endl;
            break;
        case Compression::none:
            cerr << "skipping ..." << endl;
            return 1;
        default:
            cerr << "Error." << endl;
            return 1;
    }
    return 0;
}
