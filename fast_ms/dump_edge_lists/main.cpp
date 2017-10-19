//
//  main.cpp
//  wl_node_properties
//
//  Created by denas on 10/15/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include "utils.hpp"
#include "cmd_utils.hpp"
#include "edge_list.hpp"

using namespace fdms;

void comp(InputSpec& S_fwd, const string& out_path, InputFlags& flags){
    string s = S_fwd.load_s();
    reverse_in_place(s);
    /* build stree */
    StreeOhleb<> st;
    size_type load_cst_time = load_st(st, s, S_fwd.rev_cst_fname, flags.load_stree);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;

    EdgeList e{st};
    e.write_bin(out_path);
    EdgeList g {load_edge_list_bin(out_path)};
    assert (e == g);
    
    if(flags.use_maxrep){
        sdsl::bit_vector maxrep(1);
        size_type duration = load_maxrep(maxrep, st, s, S_fwd.rev_maxrep_fname, flags.load_maxrep);
        cerr << "DONE (" << duration / 1000 << " seconds)" << endl;

        cerr << " * checking with maxrep ..." << endl;
        e.check_with_maxrep(maxrep);
        cerr << "OK" << endl;
    }
}

int main(int argc, char **argv){
    OptParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
        InputFlags flags(false, // lazy_wl
                         false,  // rank-and-fail
                         false,  // use maxrep
                         true,  // lca_parents
                         false, // space
                         false, // time
                         true,  // ans
                         false, // verbose
                         10,    // nr. progress messages for runs construction
                         10,    // nr. progress messages for ms construction
                         false, // load CST
                         false, // load MAXREP
                         1      // nthreads
                         );
        InputSpec sfwd_spec(base_dir + "mut_200s_64t_15.s");
        const string out_path = "ciao";
        comp(sfwd_spec, base_dir + out_path, flags);
    } else {
        InputFlags flags(input);
        InputSpec sfwd_spec(input.getCmdOption("-s_path"));
        const string out_path = input.getCmdOption("-out_path");
        comp(sfwd_spec, out_path, flags);
    }
    
    return 0;
}
