//
//  main.cpp
//  single_vs_double_rank
//
//  Created by denas on 1/28/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include "utils.hpp"
#include "fd_ms.hpp"
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

void double_rank_nofail(StreeOhleb<>& st, node_type v, const char_type c, const size_type ncalls){
    node_type u = v;

    for (size_type i = 0; i < ncalls; i++) {
        v = st.double_rank_nofail_wl(v, c);
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
node_type find_node(const StreeOhleb<>& st, const string s, const bool small){
    node_type v = st.root();
    size_type i = s.size() - 1;
    size_type wl_cnt = 0;

    auto block_condition = [] (node_type v, bool sm) -> bool
    {
        return (sm ? in_same_major_block(v) : !in_same_major_block(v));
    };

    i -= 1;
    v = st.single_rank_wl(v, s[i]);

    while(--i > 0){
        v = st.single_rank_wl(v, s[i]);
        wl_cnt = count_wl(st, v);
        assert (!st.is_root(v));
        assert (wl_cnt > 0);

        if(wl_cnt > 1 && wl_cnt < st.csa.sigma - 1 && block_condition(v, small))
            return v;
    }
    return st.root(); // not found
}

void dump_report_slice(StreeOhleb<>& st, string& s, size_type ntrial, size_type close, string method,
                       std::vector<size_type>& time_usage, std::vector<bool>& has_wl){
    // cout << "len_s,nwlcalls,ntrial,close,method,char,charidx,time_ms" << endl;
    for(size_type i = 0; i < time_usage.size(); i++){
        cout << s.size() << "," ;  // s_len
        cout << NWLCALLS << "," ;  // nwlcalls
        cout << ntrial << "," ;    // ntrial
        cout << close << "," ;     // is close?
        cout << method << ",";     // method name
        cout << (char) (i == 0 ? '#' : st.csa.comp2char[i]) << ","; // char
        cout << i << "," ;         // char idx
        cout << has_wl[i] << "," ;         // has wl?
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

void time_fun(StreeOhleb<>& st, string& s, const size_type ncalls, const size_type ntrials,
                      void (*fun_ptr)(StreeOhleb<>&, node_type, char_type, size_type), const string fun_name){
    node_type v = st.root();
    std::vector<size_type> run_time(st.csa.sigma);
    std::vector<bool> has_wl(st.csa.sigma);
    size_type small_interval_node[] {1, 0};


    for (auto node_closeness: small_interval_node){
        // find a node with end points in the same block
        v = find_node(st, s, node_closeness == 1);
        if(st.is_root(v))
            cerr << " ** could not find a suitable node in the ST. try a different (larger?) input" << endl;
        assert (!st.is_root(v));
        assert (count_wl(st, v) > 0);
        cerr << " ** found node (" << v.i << "," << v.j << ") with ";
        cerr << (in_same_major_block(v) ? "small" : "large") << " ";
        cerr << "(" << v.j - v.i + 1 << ") interval and ";
        cerr << count_wl(st, v) << " (of " << st.csa.sigma << " possible) Wlinks" << endl;
        for(size_type ntrial = 0; ntrial < ntrials; ntrial ++){
            cerr << " *** trial " << ntrial << " of " << ntrials << endl;
            time_method_over_chars(st, s, ncalls, v, run_time, has_wl, fun_ptr);
            dump_report_slice(st, s, ntrial + 1, node_closeness, fun_name, run_time, has_wl);
        }
    }
}




int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputFlags flags(input);

    //InputSpec s_spec("/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/datasets/testing/rnd_200_32768.s");
    //flags.load_stree = false;
    InputSpec s_spec(input.getCmdOption("-s_path"));
    string s = s_spec.load_s();
    StreeOhleb<> st;
    size_type st_time = fdms::load_st(st, s, s_spec.fwd_cst_fname, flags.load_stree);
    cerr << "DONE (" << st_time / 1000 << " seconds)" << endl;

    cout << "len_s,nwlcalls,ntrial,close,method,char,charidx,has_wl,time_ms" << endl;
    cerr << " * testing single rank" << endl;
    time_fun(st, s, NWLCALLS, (s.size() - 1 > NTRIALS ? NTRIALS : s.size() - 1),
             single_rank, "single_rank");

    cerr << " * testing double rank and fail" << endl;
    time_fun(st, s, NWLCALLS, (s.size() - 1 > NTRIALS ? NTRIALS : s.size() - 1),
             double_rank_fail, "double_rank_fail");

    cerr << " * testing double rank without fail" << endl;
    time_fun(st, s, NWLCALLS, (s.size() - 1 > NTRIALS ? NTRIALS : s.size() - 1),
             double_rank_nofail, "double_rank_no_fail");

    return 0;
}
