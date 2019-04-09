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

template<typename enc_type>
int comp(const string ms_path, const string out_suffix){
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);
    //sdsl::rrr_vector<> c_ms(ms);
    //sdsl::hyb_vector<> c_ms(ms);
    enc_type c_ms(ms);
    sdsl::store_to_file(c_ms, ms_path + out_suffix);
    return 0;
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
int comp1(const string ms_path, const string suffix){
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);

    enc_type encoder(32);
    histo_t counter;
    size_type n_runs = fill_encoder<enc_type>(ms, encoder, counter);
    vec_type c_ms(encoder, ms.size());
    std::ofstream out{ms_path + suffix, std::ios::binary};
    c_ms.writeTo(out);

    (cerr << n_runs << " runs over "
          << c_ms.getSize() << " elements ("
          << c_ms.getSize() / static_cast<float>(n_runs) << " elements / run)"
          << endl);
    for(auto item : counter)
        cout << item.first << "," << item.second << endl;
    return 0;
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
    switch(flags.compression)
    {
        case Compression::rrr:
            return comp<sdsl::rrr_vector<>>(ms_path, ms_compression::to_str(flags.compression));
        case Compression::rle:
            return comp1<CSA::RLEVector, CSA::RLEEncoder>(ms_path, ms_compression::to_str(flags.compression));
        case Compression::delta:
            return comp1<CSA::DeltaVector, CSA::DeltaEncoder>(ms_path, ms_compression::to_str(flags.compression));
        case Compression::succint:
            return comp1<CSA::SuccinctVector, CSA::SuccinctEncoder>(ms_path, ms_compression::to_str(flags.compression));
        case Compression::nibble:
            return comp1<CSA::NibbleVector, CSA::NibbleEncoder>(ms_path, ms_compression::to_str(flags.compression));
        case Compression::none:
            cerr << "skipping ..." << endl;
            break;
        default:
            cerr << "Error." << endl;
            return 1;
    }
}

