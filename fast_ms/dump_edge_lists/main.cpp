//
//  main.cpp
//  wl_node_properties
//
//  Created by denas on 10/15/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>


#include "input_spec.hpp"
#include "opt_parser.hpp"

#include "stree_sct3.hpp"
#include "cst_iterator.hpp"
#include "maxrep_vector.hpp"
#include "edge_list.hpp"


using namespace std;
using namespace fdms;


class InputFlags{
public:
    bool check, load_cst;

    InputFlags(){}
    
    InputFlags(const bool check, const bool load_cst) : check{check}, load_cst{load_cst} {}
    
    InputFlags(const InputFlags& i) : check{i.check}, load_cst{i.load_cst} {}
    
    InputFlags(const OptParser args){
        check = (args.getCmdOption("-check") == "1");
        load_cst = (args.getCmdOption("-load_cst") == "1");
    }

    void help(const string exec_name) const {
        cerr << exec_name << endl << "\t";
        cerr << "[-use_maxrep 0/1]" << endl << "\t";
        cerr << "-s_path <s path>" << endl << endl;
    }
};



void comp(InputSpec& S_fwd, InputFlags& flags){
    string s = S_fwd.load_s();
    InputSpec::reverse_in_place(s);
    
    /* build stree */
    StreeOhleb<> st;
    size_type load_cst_time = load_or_build(st, s, S_fwd.rev_cst_fname, flags.load_cst);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;

    EdgeList e{st};
    if(flags.check){
        Maxrep maxrep(st, true);
        e.check_with_maxrep(maxrep);
        cerr << "OK" << endl;
    }
    e.write_bin(S_fwd.rev_elst_fname);
}

int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec sfwd_spec;
    InputFlags flags;

    if(argc == 1){
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
        flags = InputFlags(true, false);
        sfwd_spec = InputSpec(base_dir + "mut_200s_64t_15.s");
    } else {
        flags = InputFlags(input);
        sfwd_spec = InputSpec(input.getCmdOption("-s_path"));
    }
    comp(sfwd_spec, flags);
    return 0;
}
