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

#include "utils.hpp"

using namespace fdms;



void dump(const StreeOhleb<>& st, const string fname){
    auto start = timer::now();
    cerr << " * dumping the CST to" << fname << " ";
    sdsl::store_to_file(st, fname);
    auto stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;
}

void comp(const InputSpec& s_spec, const InputFlags& flags){
    cerr << " * loading the string " << s_spec.s_fname << " ";
    auto start = timer::now();
    string s = s_spec.load_s();
    auto stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;

    StreeOhleb<> st;
    size_type load_cst_time = load_st(st, s, s_spec.fwd_cst_fname, false);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;
    dump(st, s_spec.fwd_cst_fname);


    cerr << " * reversing the string of length " << s.size() << " ";
    start = timer::now();
    reverse_in_place(s);
    stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;

    load_cst_time = load_st(st, s, s_spec.rev_cst_fname, false);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;
    dump(st, s_spec.rev_cst_fname);
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/datasets/testing/"};
        InputFlags flags(false, // lazy_wl
                         false, // rank and fail
                         false, // maxrep
                         false, // lca-parents
                         true,  // space
                         false, // time
                         true,  // ans
                         false, // verbose
                         10,    // nr. progress messages for runs construction
                         10,    // nr. progress messages for ms construction
                         false, // load CST
                         false, // load MAXREP
                         1      // nthreads
                         );
        InputSpec s_spec(base_dir + "rnd_200_64.s");
        comp(s_spec, flags);
    } else {
        InputFlags flags(input);
        InputSpec s_spec(input.getCmdOption("-s_path"));
        comp(s_spec, flags);
    }
    return 0;
}
