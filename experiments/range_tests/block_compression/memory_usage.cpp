/*
 compress a ms vector
*/

#include <iostream>
#include <fstream>
#include <vector>

//#include "sdsl/memory_management.hpp"
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


class InputFlags {
public:
    Compression compression;
    size_type start, len;

    InputFlags() { }

    InputFlags(const InputFlags& f) : compression{f.compression}{}

    InputFlags(const Compression compression) : compression{compression}{}

    InputFlags(OptParser input) {
        compression = ms_compression::parse_compression(
            input.getCmdOption("-compression")
        );
    }
};

template<typename ms_type, typename ms_sel_1_type>
void comp(const string ms_path, const InputFlags& flags){
    ms_type ms;
    sdsl::load_from_file(ms, ms_path);
    ms_sel_1_type ms_sel(&ms);

    size_type a = 0;
    for(int i=0; i<100; i++){
        a = ms_sel(1);
        a += 1;
    }
}


template<typename vec_type, typename it_type>
void comp1(const string ms_path, const InputFlags& flags){
    std::ifstream in{ms_path, std::ios::binary};
    vec_type ms(in);
    it_type* it = new it_type(ms);

    size_type a = 0;
    for(int i=0; i<100; i++){
        a = it->select(1);
        a += 1;
    }
}


int main(int argc, char **argv){
    if(argc == 1){
        (cerr << "Compress a ms vector. Creates files <ms_path>.rrr\n"
              << "Args:\n"
              << help__ms_path
              << "\t-compression <compression type>: One of: rrr, rle, delta, succint, nibble.\n"
              << endl);
        exit(0);
    }

    OptParser input(argc, argv);
    InputFlags flags(input);
    string ms_path = input.getCmdOption("-ms_path");


    switch(flags.compression)
    {
        case Compression::none:
            comp<sdsl::bit_vector, sdsl::bit_vector::select_1_type>(input.getCmdOption("-ms_path"), flags);
            break;
        case Compression::rrr:
            comp<sdsl::rrr_vector<>, sdsl::rrr_vector<>::select_1_type>(input.getCmdOption("-ms_path"), flags);
            break;
        case Compression::rle:
            comp1<CSA::RLEVector, CSA::RLEVector::Iterator>(input.getCmdOption("-ms_path"), flags);
            break;
        case Compression::delta:
            comp1<CSA::DeltaVector, CSA::DeltaVector::Iterator>(input.getCmdOption("-ms_path"), flags);
            break;
        case Compression::nibble:
            comp1<CSA::NibbleVector, CSA::NibbleVector::Iterator>(input.getCmdOption("-ms_path"), flags);
            break;
        case Compression::succint:
            comp1<CSA::SuccinctVector, CSA::SuccinctVector::Iterator>(input.getCmdOption("-ms_path"), flags);
            break;
        default:
            cerr << "Error." << endl;
            return 1;
    }
}
