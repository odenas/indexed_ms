//
//  main.cpp
//  bp
//
//  Created by denas on 11/14/16.
//  Copyright © 2016 denas. All rights reserved.
//

#include <iostream>
#include <string>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

void compute_bp(string& T){
    cst_sada<> st_of_s;
    construct_im(st_of_s, T, 1);
    cout << st_of_s.bp << endl;

    /*
    char *aa {new char[st_of_s.bp_support.size() + 1]};
    for(int i=0 ; i < st_of_s.bp_support.size(); i++)
        aa[i] = (48 + (int)st_of_s.bp[i]);
    aa[st_of_s.bp_support.size()] = '\0';
    cout << aa << endl;
    */
}

int main(int argc, const char * argv[]) {
    string T;

    if(argc == 2){ // process file
        std::ifstream in_file {argv[1]};
        if(!in_file){
            cout << "could not open file " << argv[1] << endl;
            return 1;
        }

        while (in_file >> T){
            compute_bp(T);
        }
    } else {
        string T {"abb"};
        compute_bp(T);
    }
    return 0;

}
