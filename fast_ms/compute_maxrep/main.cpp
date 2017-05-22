//
//  main.cpp
//  compute_maxrep
//
//  Created by denas on 5/9/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include <iostream>
#include <string>

#include <sdsl/suffix_trees.hpp>

#include "utils.hpp"
#include "stree_sct3.hpp"
#include "fd_ms.hpp"


using namespace std;
using namespace fdms;


typedef typename StreeOhleb<>::node_type node_type;


StreeOhleb<> st;
string s;
sdsl::bit_vector maxrep(1);


void comp(InputSpec& S_fwd, const InputFlags& flags){
    auto start = timer::now();
    cerr << " * loading and reversing the string " << S_fwd.s_fname << " ";
    s = S_fwd.load_s();
    reverse_in_place(s);
    auto stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;

    size_type load_cst_time = load_st(st, s, S_fwd.rev_cst_fname, flags.load_stree);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;

    size_type build_maxrep_time = load_maxrep(maxrep, st, s, S_fwd.rev_maxrep_fname, false);
    cerr << "DONE ( " << build_maxrep_time / 1000 << " seconds)" << endl;

    if(flags.answer){
        for(size_type i = 0; i < maxrep.size(); i++)
            cout << maxrep[i] << " ";
        cout << endl;
    }
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/datasets/testing/"};
        InputFlags flags(false, // lazy_wl
                         false, // sada cst
                         false, // maxrep
                         true, // space
                         false, // time
                         true,  // ans
                         false, // verbose
                         10,    // nr. progress messages for runs construction
                         10,    // nr. progress messages for ms construction
                         false, // load CST
                         false, // load MAXREP
                         1      // nthreads
                         );
        InputSpec tspec(base_dir + "rnd_200_64.t");
        InputSpec sfwd_spec(base_dir + "rnd_200_64.s");
        const string out_path = "0";
        comp(sfwd_spec, flags);
    } else {
        InputFlags flags(input);
        InputSpec sfwd_spec(input.getCmdOption("-s_path"));
        comp(sfwd_spec, flags);
    }
    return 0;
}
