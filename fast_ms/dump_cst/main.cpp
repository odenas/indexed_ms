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
#include <fstream>
#include <vector>

#include "input_spec.hpp"
#include "opt_parser.hpp"
#include "stree_sct3.hpp"


using namespace std;
using namespace fdms;
using timer = std::chrono::high_resolution_clock;


void dump(const StreeOhleb<>& st, const string fname){
    auto start = timer::now();
    cerr << " * dumping the CST to " << fname << " ";
    sdsl::store_to_file(st, fname);
    auto stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;
}

void comp(const InputSpec& s_spec){
    cerr << " * loading the string " << s_spec.s_fname << " ";
    auto start = timer::now();
    string s = s_spec.load_s();
    auto stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;

    StreeOhleb<> st;
    StreeOhleb<>::size_type load_cst_time = load_or_build(st, s, s_spec.fwd_cst_fname, false);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;
    dump(st, s_spec.fwd_cst_fname);

    cerr << " * reversing the string of length " << s.size() << " ";
    start = timer::now();
    InputSpec::reverse_in_place(s);
    stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;

    load_cst_time = load_or_build(st, s, s_spec.rev_cst_fname, false);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;
    dump(st, s_spec.rev_cst_fname);
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec s_spec;

    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/datasets/testing/"};
        s_spec = InputSpec(base_dir + "rnd_200_64.s");
    } else {
        s_spec = InputSpec(input.getCmdOption("-s_path"));
    }
    comp(s_spec);
    return 0;
}
