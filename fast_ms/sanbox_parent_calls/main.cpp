//
//  main.cpp
//  sanbox_parent_calls
//
//  Created by denas on 6/8/17.
//  Copyright © 2017 denas. All rights reserved.
//


#include <iostream>
#include "utils.hpp"

#include "maxrep_construction.hpp"
#include "stree_sct3.hpp"


using namespace std;
using namespace fdms;

typedef typename StreeOhleb<>::node_type node_type;
typedef void (*call_method_t) (const StreeOhleb<>& st, const node_type start_node, const char_type c, const size_type ncalls, const size_type seq_len);

size_type ncalls = 1000000;
size_type ntrials = 5;
size_type max_seq_len = 5;
vector<std::pair<node_type, char_type>> starting_points(max_seq_len);



void call_pseq(const StreeOhleb<>& st, const node_type start_node, const char_type c, const size_type ncalls, const size_type seq_len){
    node_type v = start_node;
    //cout << v.i << "," << v.j << endl;
    for(size_type i = 0; i < ncalls; i++){
        for(size_type j = 0; j < seq_len; j++)
            v = st.parent_sequence(v, c);
    }
    //cout << v.i << "," << v.j << endl;
}

void call_lca(const StreeOhleb<>& st, const node_type start_node, const char_type c, const size_type ncalls, const size_type seq_len){
    node_type v = start_node;
    //cout << v.i << "," << v.j << endl;
    for(size_type i = 0; i < ncalls; i++)
        v = st.maxrep_ancestor(v, c);
    //cout << v.i << "," << v.j << endl;
}



//bool in_same_major_block(const node_type v) { return ((v.i>>8) == (v.j>>8)); }


bool has_wl(const StreeOhleb<>& st, const node_type v, const char_type c){
    return !st.is_root(st.double_rank_fail_wl(v, c));
}

size_type parent_depth(const StreeOhleb<>& st, const node_type start_node, const char_type c){
    size_type d = 0;
    node_type u = start_node, v = st.double_rank_fail_wl(u, c);

    while((!st.is_root(u)) and st.is_root(v)){
        u = st.parent(u);
        v = st.double_rank_fail_wl(u, c);
        d += 1;
    }
    return d;
}

std::pair<node_type, char_type> starting_point(const StreeOhleb<>& st, const string& s, const size_type seq_len){
    node_type v = st.root();
    size_type k = s.size();

    // do 2*seq_len wl() calls
    for (size_type i=0; i < 2*seq_len; i++, --k) {
        v = st.double_rank_fail_wl(v, s[k]);
    }

    while(--k > 1){
        // check if valid for some symbol
        for(size_type i = 1; i < st.csa.sigma; i++){
            char_type c = st.csa.comp2char[i];
            if(parent_depth(st, v, c) == seq_len)
                return std::make_pair(v, c);
        }
        v = st.double_rank_fail_wl(v, s[--k]);
    }
    cerr << "could not find node with seq_len " << seq_len << endl;
    return std::make_pair(st.root(), st.csa.comp2char[0]);
}

std::pair<size_type, size_type> i_width(const StreeOhleb<>& st, const node_type starting_node, const size_type seq_len){
    node_type u = starting_node;
    size_type width = u.j - u.i + 1;

    for(size_type i = 0; i < seq_len; i++)
        u = st.parent(u);
    return std::make_pair(width,  u.j - u.i + 1);
}

void time_fun(StreeOhleb<>& st, string& s, const size_type ncalls, const size_type ntrials,
              call_method_t call_ptr, const string fun_name, const size_type max_seq_len,
              const vector<std::pair<node_type, char_type>>& starting_points){
    for(size_type seq_len = 1; seq_len <= max_seq_len; seq_len++){
        node_type v = starting_points[seq_len].first;
        char_type c = starting_points[seq_len].second;
        std::pair<size_type, size_type> width = i_width(st, v, seq_len);

        for(size_type i = 0; i < ntrials; i++){
            auto start_time = timer::now();
            call_ptr(st, v, c, ncalls, seq_len);
            auto end_time = timer::now();
            size_type run_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

            // cout << "ntrial,cnt,seq_len,method,char,char_idx,start_width,end_width,value"
            cout << i + 1<< "," ;                    // ntrial
            cout << ncalls << "," ;               // ncalls
            cout << seq_len << "," ;              // seq_len
            cout << fun_name << "," ;             // method
            cout << c << "," ;                    // char
            cout << (int)st.csa.char2comp[c] << "," ;  // char_idx
            cout << width.first << ",";             // width of starting node
            cout << width.second << ",";             // width of ending node
            cout << run_time << endl;             // value
        }
    }
}


int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputFlags flags(input);

    //InputSpec s_spec("/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/datasets/testing/rnd_200_256.s");
    //flags.load_stree = false;
    InputSpec s_spec(input.getCmdOption("-s_path"));
    string s = s_spec.load_s();
    StreeOhleb<> st;
    size_type st_time = fdms::load_st(st, s, s_spec.fwd_cst_fname, flags.load_stree);
    cerr << "DONE (" << st_time / 1000 << " seconds)" << endl;

    starting_points[0] = std::make_pair(st.root(), st.csa.comp2char[0]);
    for (size_type i = 1; i <= max_seq_len; i++) {
        std::pair<node_type, char_type> start_at = starting_point(st, s, i);
        if(st.is_root(start_at.first)){
            cerr << " ** [" << i << "] coult nod find node. (" << start_at.first.i << "," << start_at.first.j << "), char = " << start_at.second << endl;
            exit(1);
        }
        assert (!st.is_root(start_at.first));
        cerr << " ** [" << i << "] found node (" << start_at.first.i << "," << start_at.first.j << "), char = " << start_at.second << endl;
        starting_points[i] = start_at;
    }

    vector<std::pair<call_method_t, string>> methods = {std::make_pair(call_pseq, "pseq"), std::make_pair(call_lca, "lca")};
    cout << "ntrial,cnt,seq_len,method,char,char_idx,start_width,end_width,value" << endl;
    for(auto item: methods){
        cerr << " * testing " << item.second << endl;
        time_fun(st, s, ncalls, ntrials, item.first, item.second, max_seq_len, starting_points);
    }
    return 0;
}
