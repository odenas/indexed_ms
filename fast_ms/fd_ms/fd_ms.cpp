/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */

#include <iostream>
#include <string>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bp_support.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/bit_vectors.hpp>

#include "fd_ms.hpp"


using namespace std;
//using namespace sdsl;
using namespace fdms;

Mstat compute_ms(string& T,
                 string& Sfwd, string& Sfwdbp,
                 string& Srev, string& Srevbp,
                 const bool verbose){

    bit_vector bfwd(Sfwdbp.size()), brev(Srevbp.size());
    for(int i=0; i<Sfwdbp.size(); i++){
        bfwd[i] = ((unsigned char)Sfwdbp[i] - 48);
        brev[i] = ((unsigned char)Srevbp[i] - 48);
    }
    fdms::bp_support_sada<> Bpsfwd(&bfwd);
    fdms::bp_support_sada<> Bpsrev(&brev);

    if(verbose && false){
        cst_sct3<> st_of_s, st_of_s_rev;
        construct_im(st_of_s, Sfwd, 1);
        construct_im(st_of_s_rev, Srev, 1);

        cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
        csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s.csa);
        cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
        csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s_rev.csa);
    }

    Bwt Bwtfwd((const unsigned char *)Sfwd.c_str());
    for(size_type i=1; i <= Bwtfwd.bwt_len; i++)
        cout << "rank(a, " << i << ") = " << (int) Bwtfwd.rank(i, 'a') << endl;
    cout << endl;
    
    Bwt Bwtrev((const unsigned char *)Srev.c_str());

    Mstat MS(T,
             Sfwd, Bwtfwd, Bpsfwd,
             Srev, Bwtrev, Bpsrev,
             verbose);

    return MS;
}


int main(int argc, char **argv){
    string T, S, Sbp, Srev, Srevbp;

    if(argc == 2){ // process file
        std::ifstream in_file {argv[1]};
        if(!in_file){
            cout << "could not open file " << argv[1] << endl;
            return 1;
        }

        while (in_file >> T >> S >> Sbp >> Srev >> Srevbp){
            cout << T + " " + S + " ";
            compute_ms(T, S, Sbp, Srev, Srevbp, 0).dump();
        }
    } else {
        string T {"aaabbba"};
        string S {"aabbaba"};
        string Sbp {"11011010110100011101001000"};
        string Srev {"ababbaa"};
        string Srevbp {"11011010110100011101001000"};

        memory_monitor::start();
        compute_ms(T, S, Sbp, Srev, Srevbp, 1).dump();
        memory_monitor::stop();
        cout << memory_monitor::peak() << " bytes." << endl;

        //cout << "usage: " << argv[0] << "<file of strings>" << endl;
        //cout << "or" << endl;
        //cout << "usage: " << argv[0] << "<T> <S>" << endl;
        //return 1;
    }

    return 0;
}

