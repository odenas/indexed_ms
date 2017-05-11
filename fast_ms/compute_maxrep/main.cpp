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


void comp(InputSpec& T, InputSpec& S_fwd, const string& out_path, InputFlags& flags){
    auto start = timer::now();
    cerr << " * loading and reversing the string " << S_fwd.s_fname << " ";
    s = S_fwd.load_s();
    reverse_in_place(s);
    auto stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;

    start = timer::now();
    if(flags.load_stree){
        cerr << " * loading the CST T(s') from " << S_fwd.rev_cst_fname << " ";
        sdsl::load_from_file(st, S_fwd.rev_cst_fname);
    } else {
        cerr << " * building the CST T(s') of length " << s.size() << " ";
        sdsl::construct_im(st, s, 1);
    }
    stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;

    maxrep.resize(s.size() + 1); sdsl::util::set_to_value(maxrep, 0);
    start = timer::now();
    cerr << " * computing MAXREP ";
    build_maxrep_ohleb(st, maxrep);
    stop = timer::now();
    cerr << "DONE ( " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " milliseconds)" << endl;

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
                         1      // nthreads
                         );
        InputSpec tspec(base_dir + "rnd_200_64.t");
        InputSpec sfwd_spec(base_dir + "rnd_200_64.s");
        const string out_path = "0";
        comp(tspec, sfwd_spec, out_path, flags);
    } else {
        InputFlags flags(input);
        InputSpec tspec(input.getCmdOption("-t_path"));
        InputSpec sfwd_spec(input.getCmdOption("-s_path"));
        const string out_path = input.getCmdOption("-out_path");
        comp(tspec, sfwd_spec, out_path, flags);
    }
    return 0;
}
