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
#include "fd_ms/p_ms_vector.hpp"
#include "fd_ms/stree_sct3.hpp"

#include "rlcsa/bits/bitbuffer.h"

#include "../malloc_count/malloc_count.h"

using namespace fdms;
using namespace std;

typedef StreeOhleb<> cst_t;
typedef typename cst_t::size_type size_type;
typedef typename ms_compression::compression_types Compression;
typedef typename sdsl::int_vector<64> vec_type;

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
        if (compression != Compression::delta and compression != Compression::rle){
            throw string{"only rle or delta compression supported. "};
        }
    }
};


void write_through(CSA::WriteBuffer& buff, std::ofstream& out, const size_type value){
    if(!buff.writeDeltaCode(value)){
        buff.writeTo(out);
        buff.reset();
        buff.writeDeltaCode(value);
    }
}

size_type delta_encode(const string ms_path){
    vec_type freq;
    sdsl::load_from_file(freq, ms_path);
    CSA::WriteBuffer buff(100);

    // writing just delta encoding without RLE
    std::ofstream out{ms_path + ".delta", std::ios::binary};
    for(size_type i = 0; i < 30; i++)
        write_through(buff, out, freq[i]);
    buff.writeTo(out);
}

size_type rle_encode(const string ms_path){
    vec_type freq;
    sdsl::load_from_file(freq, ms_path);
    CSA::WriteBuffer buff(100);

    std::ofstream out{ms_path + ".rle", std::ios::binary};
    size_type previous_frequency = freq[0];
    size_type run_length = 0;

    for(size_type i = 0; i < 30; i++){
        if(freq[i] == previous_frequency)
            run_length += 1;
        else{
            write_through(buff, out, previous_frequency);
            write_through(buff, out, run_length);
            previous_frequency = freq[i];
            run_length = 1;
        }
    }
    write_through(buff, out, previous_frequency);
    write_through(buff, out, run_length);
    buff.writeTo(out);
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputFlags flags;

    if(argc == 1){
        (cerr << "Compress a freq vector with rle-encoding. Creates files <path>.rle\n"
              << "Args:\n"
              << help__freq_path
              << "\t-compression <compression type>: One of: rle, delta." << endl
              << endl);
        exit(0);
    }
    string freq_path = input.getCmdOption("-freq_path");
    try{
        flags = InputFlags(input);
        switch(flags.compression)
        {
            case Compression::rle:
                cout << rle_encode(freq_path) << endl;
                break;
            case Compression::delta:
                cout << delta_encode(freq_path) << endl;
                break;
            default:
                cerr << "Error." << endl;
                return 1;
        }
    } catch(string s){
        throw s;
    }
    return 0;
}
