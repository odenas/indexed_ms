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
    size_t sample_freq;

    InputFlags(){}
    
    InputFlags(const bool check, const bool load_cst, const size_t sample_freq) :
        check{check}, load_cst{load_cst}, sample_freq{sample_freq} {}

    InputFlags(const InputFlags& i) : check{i.check}, load_cst{i.load_cst}, sample_freq{i.sample_freq} {}
    
    InputFlags(const OptParser args){
        check = (args.getCmdOption("-check") == "1");
        load_cst = (args.getCmdOption("-load_cst") == "1");
        sample_freq = (static_cast<size_t>(std::stoi(args.getCmdOption("-sample_freq"))));
    }
};



void comp(InputSpec& S_fwd, InputFlags& flags){
    string s = S_fwd.load_s();
    InputSpec::reverse_in_place(s);
    
    /* build stree */
    StreeOhleb<> st;
    size_type load_cst_time = load_or_build(st, s, S_fwd.rev_cst_fname, flags.load_cst);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;

    EdgeList e{st, flags.sample_freq};
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
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/code_test/maxrep_inputs/"};
        flags = InputFlags(true, false, 2);
        sfwd_spec = InputSpec(base_dir + "rnd_20s_dis_10t_abcd.s");
    } else {
        flags = InputFlags(input);
        sfwd_spec = InputSpec(input.getCmdOption("-s_path"));
    }
    comp(sfwd_spec, flags);
    return 0;
}
