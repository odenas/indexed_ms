//
//  main.cpp
//  sanbox_parent_calls
//
//  Created by denas on 6/8/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/parent_depth_list.hpp"




using namespace std;
using namespace fdms;


typedef typename StreeOhleb<>::node_type node_type;
typedef typename StreeOhleb<>::char_type char_type;
typedef void (*call_method_t) (const StreeOhleb<>& st, const node_type start_node, const char_type c, const size_type seq_len);



class InputFlags{
public:
    bool load_cst;
    size_t repeat;
    
    InputFlags(){}
    
    InputFlags(const bool load_cst, size_t repeat) : load_cst{load_cst}, repeat{repeat} {}
    
    InputFlags(const InputFlags& i) : load_cst{i.load_cst}, repeat{i.repeat} {}
    
    InputFlags(const OptParser args){
        load_cst = (args.getCmdOption("-load_cst") == "1");
        repeat = (static_cast<size_t>(std::stoi(args.getCmdOption("-repeat"))));
    }
};


void call_pseq(const StreeOhleb<>& st, const node_type start_node, const char_type c, const size_type seq_len){
    node_type v = start_node;
    //cout << v.i << "," << v.j << endl;
    for(size_type j = 0; j < seq_len; j++)
        v = st.parent_sequence(v, c);
    //cout << v.i << "," << v.j << endl;
}


void call_lca(const StreeOhleb<>& st, const node_type start_node, const char_type c, const size_type seq_len){
    node_type v = start_node;
    //cout << v.i << "," << v.j << endl;
    v = st.maxrep_ancestor(v, c);
    //cout << v.i << "," << v.j << endl;
}

vector<node_with_depth> constant_depth_nodes(vector<node_with_depth> u, const size_type depth){
    vector<node_with_depth> v;
    for(auto nwd : u){
        if(nwd.m_depth == depth)
            v.push_back(nwd);
    }
    return v;
}


size_type time_fun(const StreeOhleb<>& st, const vector<node_with_depth> v, call_method_t call_ptr, const size_type ntrials){
    auto start_time = timer::now();
    for(size_type ntrial = 0; ntrial < ntrials; ntrial++){
        for(auto nwd : v)
            call_ptr(st, nwd.m_node, nwd.m_c, nwd.m_depth);
    }
    auto end_time = timer::now();

    return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
}

void report(const size_type depth, const string method, const size_type nwd_cnt, const size_type run_time){
    cout << depth << "," ;             // parent_depth
    cout << method << "," ;          // method
    cout << nwd_cnt << "," ;          // nwd_cnt
    cout << run_time << endl;          // value
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputSpec sfwd_spec;
    InputFlags flags;

    if(argc == 1){
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/code_test/maxrep_inputs/"};
        sfwd_spec = InputSpec(base_dir + "rnd_20s_dis_10t_abcd.s");
        flags.repeat = 10;
    } else {
        flags = InputFlags(input);
        sfwd_spec = InputSpec(input.getCmdOption("-s_path"));
    }
    string s = sfwd_spec.load_s();
    InputSpec::reverse_in_place(s);
    
    StreeOhleb<> st;
    size_type st_time = load_or_build(st, s, sfwd_spec.rev_cst_fname, flags.load_cst);
    cerr << "DONE (" << st_time / 1000 << " seconds)" << endl;
    
    NwdList nlst = load_nwd_list_bin(sfwd_spec.fwd_nwdlst_fname);
    cerr << nlst.repr() << endl;

    cout << "parent_depth,method,nwd_cnt,value" << endl;
    for(size_type parent_depth = 1; parent_depth < nlst.max_depth1; parent_depth++){
        vector<node_with_depth> v = constant_depth_nodes(nlst.m_vec1, parent_depth);
        std::random_shuffle(v.begin(), v.end());
        report(parent_depth, "pseq", v.size(), time_fun(st, v, call_pseq, flags.repeat));
        report(parent_depth, "lca", v.size(), time_fun(st, v, call_lca, flags.repeat));
    }
    return 0;
}
