//
//  main.cpp
//  sandbox_test_maxrep
//
//  Created by denas on 5/20/17.
//  Copyright © 2017 denas. All rights reserved.
//


#include <iostream>
#include "utils.hpp"

#include "maxrep_construction.hpp"
#include "stree_sct3.hpp"


#define NWLCALLS 1000000
#define NTRIALS  5

using namespace std;
using namespace fdms;

typedef typename StreeOhleb<>::node_type node_type;


void single_rank(StreeOhleb<>& st, node_type v, const char_type c, const size_type ncalls){
    node_type u = v;

    for (size_type i = 0; i < ncalls; i++) {
        v = st.single_rank_wl(v, c);
        v = u;
    }
}

void double_rank_fail(StreeOhleb<>& st, node_type v, const char_type c, const size_type ncalls){
    node_type u = v;

    for (size_type i = 0; i < ncalls; i++) {
        v = st.double_rank_fail_wl(v, c);
        v = u;
    }
}

void double_rank_maxrep(StreeOhleb<>& st, node_type v, const char_type c, const size_type ncalls){
    node_type u = v;

    for (size_type i = 0; i < ncalls; i++) {
        v = st.double_rank_nofail_wl(v, c);
        v = u;
    }
}


size_type count_wl(const StreeOhleb<>& st, const node_type v){
    size_type wl_count = 0;
    for(size_type cidx = 0; cidx < st.csa.sigma; cidx++){
        if(st.is_root(st.single_rank_wl(v, st.csa.comp2char[cidx])))
            wl_count += 1;
    }
    return wl_count;
}

bool in_same_major_block(const node_type v) { return ((v.i>>8) == (v.j>>8)); }

/* small = in same major block */
node_type find_node(const StreeOhleb<>& st, const sdsl::bit_vector maxrep, const string s, const bool small, const bool maximal){
    auto block_condition = [] (node_type v, bool sm) -> bool
    {
        return (sm ? in_same_major_block(v) : !in_same_major_block(v));
    };

    auto wlcount_condition = [&] (const node_type v) -> bool
    {
        size_type wl_cnt = count_wl(st, v);
        assert (!st.is_root(v));
        return (wl_cnt > 0 && wl_cnt < st.csa.sigma);
    };

    auto maximal_condition = [&] (const node_type v) -> bool
    {
        bool is_maximal = (maxrep[v.i] == 1 && maxrep[v.j] == 1);
        return (maximal ? is_maximal : !is_maximal);
    };


    // navigate the tree
    node_type currnode = st.root(), nextnode = st.root();
    bool direction_down = true;

    do{
        if(direction_down){
            if((!st.is_leaf(currnode)) && (!st.is_root(currnode))){ // process currnode
                if(wlcount_condition(currnode) && maximal_condition(currnode) && block_condition(currnode, small))
                    return currnode;
            }

            nextnode = st.first_child(currnode);
            if(st.is_root(nextnode))
                direction_down = false;
            else
                currnode = nextnode;
        } else {
            nextnode = st.sibling(currnode);
            if(st.is_root(nextnode)) {
                currnode = st.parent(currnode);
            } else {
                currnode = nextnode;
                direction_down = true;
            }
        }
    } while(!st.is_root(currnode));

    return st.root(); // not found
}

void dump_report_slice(StreeOhleb<>& st, string& s, size_type ntrial, size_type close, string method,
                       std::vector<size_type>& time_usage, std::vector<bool>& has_wl, const size_type is_max){
    // cout << "len_s,nwlcalls,ntrial,close,method,char,charidx,time_ms" << endl;
    for(size_type i = 0; i < time_usage.size(); i++){
        cout << s.size() << "," ;  // s_len
        cout << NWLCALLS << "," ;  // nwlcalls
        cout << ntrial << "," ;    // ntrial
        cout << close << "," ;     // is close?
        cout << method << ",";     // method name
        cout << (char) (i == 0 ? '#' : st.csa.comp2char[i]) << ","; // char
        cout << i << "," ;         // char idx
        cout << has_wl[i] << "," ; // has wl?
        cout << is_max << ",";     // is_maximal
        cout << time_usage[i] << endl;
    }
}


void time_method_over_chars(StreeOhleb<>& st, string& s, const size_type ncalls, const node_type v,
                            std::vector<size_type>& run_time, std::vector<bool>& has_wl,
                            void (*fun_ptr)(StreeOhleb<>&, node_type, char_type, size_type)){
    for(size_type cidx = 0; cidx < st.csa.sigma; cidx++){
        char_type c = st.csa.comp2char[cidx];
        has_wl[cidx] = !(st.is_root(st.single_rank_wl(v, c)));

        auto start_time = timer::now();
        fun_ptr(st, v, c, ncalls);
        run_time[cidx] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
    }

}

void time_fun(StreeOhleb<>& st, string& s, sdsl::bit_vector& maxrep, const size_type ncalls, const size_type ntrials,
              void (*fun_ptr)(StreeOhleb<>&, node_type, char_type, size_type), const string fun_name){
    node_type v = st.root();
    std::vector<size_type> run_time(st.csa.sigma);
    std::vector<bool> has_wl(st.csa.sigma);
    size_type small_interval_node[] {1, 0};
    size_type is_maxrepeat_node[] {1, 0};



    for(auto maximality: is_maxrepeat_node){
        for (auto node_closeness: small_interval_node){
            // find a node with end points in the same block
            v = find_node(st, maxrep, s, node_closeness == 1, maximality == 1);
            if(st.is_root(v)){
                cerr << " ** could not find a suitable node in the ST";
                cerr << " ** maximality = " << maximality << ", interval = " << node_closeness;
                cerr << ". try a different (larger?) input" << endl;
            }
            assert (!st.is_root(v));
            assert (count_wl(st, v) > 0);

            cerr << " ** found node (" << v.i << "," << v.j << "), ";
            cerr << (maxrep[v.i] == 1 && maxrep[v.j] == 1 ? "maximal" : "non-maximal") << ", ";
            cerr << (in_same_major_block(v) ? "small" : "large") << " ";
            cerr << "(" << v.j - v.i + 1 << ") interval and ";
            cerr << count_wl(st, v) << " (of " << st.csa.sigma << " possible) Wlinks" << endl;

            for(size_type ntrial = 0; ntrial < ntrials; ntrial ++){
                cerr << " *** trial " << ntrial << " of " << ntrials << endl;
                time_method_over_chars(st, s, ncalls, v, run_time, has_wl, fun_ptr);
                dump_report_slice(st, s, ntrial + 1, node_closeness, fun_name, run_time, has_wl, maximality);
            }
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


    sdsl::bit_vector maxrep(s.size() + 1);
    size_type mrep_time = fdms::load_maxrep(maxrep, st, s, s_spec.rev_maxrep_fname, flags.load_maxrep);
    cerr << "DONE ( " << mrep_time << " milliseconds)" << endl;

    cout << "len_s,nwlcalls,ntrial,close,method,char,charidx,has_wl,is_maximal,time_ms" << endl;
    cerr << " * testing single rank" << endl;
    time_fun(st, s, maxrep, NWLCALLS, (s.size() - 1 > NTRIALS ? NTRIALS : s.size() - 1),
             single_rank, "single_rank");

    cerr << " * testing double rank and fail" << endl;
    time_fun(st, s, maxrep, NWLCALLS, (s.size() - 1 > NTRIALS ? NTRIALS : s.size() - 1),
             double_rank_fail, "double_rank_fail");

    cerr << " * testing double rank with maxrep" << endl;
    time_fun(st, s, maxrep, NWLCALLS, (s.size() - 1 > NTRIALS ? NTRIALS : s.size() - 1),
             double_rank_maxrep, "double_rank_maxrep");
    
    return 0;
}
