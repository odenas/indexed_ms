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


size_type compute_length(const size_type n_ones, const size_type n_zeros){
    if (n_ones < n_zeros)
        return 2 * (n_zeros - n_ones);
    else if (n_ones < 2 * n_zeros)
        return 2 * (n_ones - n_zeros) + 1;
    return n_ones;
}

template<typename enc_type>
size_type fill_encoder_runCompression(sdsl::bit_vector ms, enc_type& encoder) {
    size_type nBitsInOutput, nOnes, nZeros, length, nRuns;
    
    nBitsInOutput = 0; nZeros = 0; nOnes = 0; nRuns = 0;
    for (size_type i=0; i<ms.size(); i++) {
        if (ms[i] == 1) {
            nOnes++;
            continue;
        }

        // ms[i] = 0
        if (nOnes == 0) {
            nZeros++;
            continue;
        }

        // ms[i] == 0 following a sequence of nOnes 1s.
        nRuns += 1;
        length = compute_length(nOnes, nZeros);

        encoder.addRun(nBitsInOutput,length);
        nOnes = 0;
        nZeros = 1;
        nBitsInOutput += (length + 1);
    }
    if (nOnes > 0) {
        nRuns += 1;
        length = compute_length(nOnes, nZeros);
        encoder.addRun(nBitsInOutput, length);
        nBitsInOutput += length;
    }
    encoder.flush();
    printf("bitsInMS=%lu bitsInCompressedMS=%lu \n",ms.size(),nBitsInOutput);
    return nRuns;
}





template<typename vec_type, typename enc_type>
size_type comp1(const string ms_path, const InputFlags& flags){
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms, ms_path);

    size_type from = abs_point();
    enc_type encoder(32);
    histo_t counter;
    
//    size_type n_runs = fill_encoder<enc_type>(ms, encoder, counter);    
    size_type n_runs=fill_encoder_runCompression<enc_type>(ms,encoder);
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



void runsHistogram(const string ms_path) {
    size_type i, nZeros, nOnes, delta;
    unsigned int MAX_DELTA = 1000;
    size_type histogram[MAX_DELTA+1+MAX_DELTA];
    
    sdsl::bit_vector ms;
    sdsl::load_from_file(ms,ms_path);
    for (i=0; i<MAX_DELTA+1+MAX_DELTA; i++) histogram[i]=0;

    // Computing histogram
    nZeros=0; nOnes=0;
        for (i=0; i<ms.size(); i++) {
            if (ms[i]==0) {
            if (nOnes==0) nZeros++;
            else {
                delta=MAX_DELTA+(nZeros-nOnes);
                if (delta>MAX_DELTA+MAX_DELTA) delta=MAX_DELTA+MAX_DELTA;
                else if (delta<0) delta=0;
                histogram[delta]++;
                nZeros=1; nOnes=0;
            }
        }
            else nOnes++;
        }
    delta=MAX_DELTA+(nZeros-nOnes);
    if (delta>MAX_DELTA+MAX_DELTA) delta=MAX_DELTA+MAX_DELTA;
    else if (delta<0) delta=0;
    histogram[delta]++;
    
    // Printing histogram
    for (i=0; i<MAX_DELTA+1+MAX_DELTA; i++) printf("%lu,%lu \n",(-MAX_DELTA+i),histogram[i]);
}




void runsCorrelation(const string ms_path) {
    const size_type HISTORY_SIZE = 100;
        int i, j;
    size_type nZeros, nOnes, nextHistory;
    sdsl::bit_vector ms;
    size_type history[HISTORY_SIZE];
    
    sdsl::load_from_file(ms,ms_path);
    for (i=0; i<HISTORY_SIZE; i++) history[i]=0;
    nZeros=0; nOnes=0; nextHistory=0;
        for (i=0; i<ms.size(); i++) {
            if (ms[i]==0) {
            if (nOnes==0) nZeros++;
            else {
                history[nextHistory]=nZeros;
                for (j=nextHistory; j>=0; j--) printf("%lu,",history[j]);
                for (j=HISTORY_SIZE-1; j>nextHistory; j--) printf("%lu,",history[j]);
                printf("\n");
                nextHistory=(nextHistory+1)%HISTORY_SIZE;
                nZeros=1; nOnes=0;
            }
        }
            else nOnes++;
        }
    history[nextHistory]=nZeros;
    for (j=nextHistory; j>=0; j--) printf("%lu,",history[j]);
    for (j=HISTORY_SIZE-1; j>nextHistory; j--) printf("%lu,",history[j]);
    printf("\n");
}    





int main(int argc, char **argv){
    OptParser input(argc, argv);
    string ms_path;
    InputFlags flags;

    if(argc == 1){
        (cerr << "Compress a ms vector. Creates files <ms_path>.xxx\n"
              << "Args:\n"
              << help__ms_path
              << "\t-compression <compression type>: One of: rrr, hybrid, rle, delta, succint, nibble.\n"
              << endl);
        exit(0);
    }
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
        runsCorrelation(ms_path);
//runsHistogram(ms_path); 
        //cout << comp1<CSA::SuccinctVector, CSA::SuccinctEncoder>(ms_path, flags) << endl;
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
