//
//  main.cpp
//  dump_maxrep
//
//  Created by denas on 5/22/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include <iostream>
#include <string>

#include <sdsl/suffix_trees.hpp>

#include "utils.hpp"
#include "cmd_utils.hpp"
#include "stree_sct3.hpp"


using namespace std;
using namespace fdms;


typedef typename StreeOhleb<>::node_type node_type;


StreeOhleb<> st;
string s;
sdsl::bit_vector maxrep(1);


void comp(const InputSpec& s_spec, const InputFlags& flags){
    auto start = timer::now();
    cerr << " * loading and reversing the string " << s_spec.s_fname << " ";
    s = s_spec.load_s();
    reverse_in_place(s);
    auto stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;

    size_type load_cst_time = load_st(st, s, s_spec.rev_cst_fname, flags.load_stree);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;

    maxrep.resize(s.size() + 1); sdsl::util::set_to_value(maxrep, 0);
    start = timer::now();
    cerr << " * computing MAXREP ";
    build_maxrep_ohleb<StreeOhleb<>, sdsl::bit_vector, StreeOhleb<>::node_type>(st, maxrep);
    stop = timer::now();
    cerr << "DONE ( " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " milliseconds)" << endl;


    cerr << " * dumping MAXREP to " << s_spec.rev_maxrep_fname << " ";
    start = timer::now();
    sdsl::store_to_file(maxrep, s_spec.rev_maxrep_fname);
    stop = timer::now();
    cerr << "DONE ( " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " milliseconds)" << endl;
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
        InputFlags flags(false, // lazy_wl
                         true,  // rank-and-fail
                         false,  // use maxrep
                         false,  // lca_parents
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
        InputSpec s_spec(base_dir + "mut_200s_64t_15.s");
        comp(s_spec, flags);
    } else {
        InputFlags flags(input);
        InputSpec s_spec(input.getCmdOption("-s_path"));
        comp(s_spec, flags);
    }
    return 0;
}
