//
//  main.cpp
//  dump_cst
//
//  Created by denas on 1/23/17.
//  Copyright © 2017 denas. All rights reserved.
//

/*
 dump a CST to a file. to be used for later use
 */

#include <iostream>

#include "stree_sct3.hpp"
#include "utils.hpp"

using namespace fdms;

template<typename tree_type>
void dump_stree(const string s, const string out_path){
    cerr << "building T(s) of length ..." << s.size() << endl;
    tree_type st;
    sdsl::construct_im(st, s, 1);
    cerr << "dumping to  " << out_path << endl;
    sdsl::store_to_file(st, out_path);
}


int main(int argc, char  **argv) {
    if(argc == 1){
        cout << argv[0] << " -s_path" << endl;
        cout << " running on a sample string ..." << endl;

        const string s_path = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/datasets/testing/abcde200_1024s.txt"};
        fdms::InputSpec sfwd_spec(s_path);
        string s = sfwd_spec.load_s();
        fdms::StreeOhleb<> st;
        sdsl::construct_im(st, s, 1);
        sdsl::store_to_file(st, sfwd_spec.s_fname + ".fwd.stree");
    } else {
        fdms::InputSpec sfwd_spec(argv[1]);

        string s = sfwd_spec.load_s();
        dump_stree<StreeOhleb<>>(s, sfwd_spec.s_fname + ".fwd.stree");

        reverse_in_place(s);
        dump_stree<StreeOhleb<>>(s, sfwd_spec.s_fname + ".rev.stree");
    }
}
