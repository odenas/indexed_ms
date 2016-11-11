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
using namespace sdsl;

Mstat compute_ms(string& T, string& S, const bool verbose){
    string Srev {S};
    for(int i=0; i<S.length(); i++)
        Srev[S.length() - i - 1] = S[i];

    if(verbose && false){
        cst_sct3<> st_of_s, st_of_s_rev;
        construct_im(st_of_s, S, 1);
        construct_im(st_of_s_rev, Srev, 1);

        cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
        csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s.csa);
        cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
        csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s_rev.csa);
    }

    Bwt bwt_fw((const unsigned char *)S.c_str());
    Bwt bwt_rev((const unsigned char *)Srev.c_str());
    Mstat MS(T, S, bwt_fw, Srev, bwt_rev, verbose);
    return MS;
}


int main(int argc, char **argv){
    string T, S, BWTfw, BWTrev, A, Cstr;

    if(argc == 2){ // process file
        std::ifstream in_file {argv[1]};
        if(!in_file){
            cout << "could not open file " << argv[1] << endl;
            return 1;
        }

        while (in_file >> T >> S >> BWTfw >> BWTrev >> A >> Cstr){
            cout << T + " " + S + " ";
            compute_ms(T, S, 0).dump();
        }
    } else {
        string T {"bbbabbababaaabbbbabaaababbbbabaaaabaabbaaaabbaaaaaabababbbaaabbaaaabaabbababbabbabaaaaaaaabaababbababaaabaaaaababbbabbbaaaabbabbbbbaabbbbbbabaabbbbabbabbbaaabbbbbaaaaaabbbaaaabbaaabaaabaaabbbabaabbaaababaababbaaaababaabababbbaaaaabbbaabaabaababbbababbaabaa"};
        string S {"bbabbabaababbabbbaabbaabbaaaaaabaaabbbbbbbaaaaabbabaaababababbbabbaabaaaaaabbaaaabbabbbbbbbaabaaaaababbabbbabaaaaabbabbabbababababbaabbbaaabaaaaaabbbababbbaabbbaaabaaaaabbababababbaaaabbabbbbaaabaabbbaabaaabbaabbaabbbaaaabbabbabbabbaaabaababbbbababbbbabaab"};

//        std::fstream in_file {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/a.txt"};
//        in_file >> T >> S >> BWTfw >> BWTrev >> A >> Cstr;

        memory_monitor::start();
        compute_ms(T, S, 1).dump();
        memory_monitor::stop();
        cout << memory_monitor::peak() << " bytes." << endl;

        //cout << "usage: " << argv[0] << "<file of strings>" << endl;
        //cout << "or" << endl;
        //cout << "usage: " << argv[0] << "<T> <S>" << endl;
        //return 1;
    }

    return 0;
}

