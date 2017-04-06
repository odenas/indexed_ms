//
//  main.cpp
//  lazy_wl_performance
//
//  Created by denas on 1/28/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include "basic.hpp"
#include "fd_ms.hpp"
#include "stree_sct3.hpp"

using namespace std;
using namespace fdms;
using timer = std::chrono::high_resolution_clock;

typedef pair<size_type, size_type> IInterval;
typedef map<std::string, size_type> Counter;
typedef typename StreeOhleb<>::node_type node_type;

void benchmark_lazywl(string& s_rev, StreeOhleb<>& st, const size_type ntrials, size_type trial_length){
    size_type nt = 0;
    size_type i = 0;
    size_type k = s_rev.size() - 1;
    node_type v = st.root();

    while(nt++ < ntrials){
        for(i = 0; i < trial_length; i++)
            v = st.lazy_wl(v, s_rev[k--]);
        if(i > 0) // finish completing the new node
            st.lazy_wl_followup(v);

        if(k < ntrials)
            k = s_rev.size() - 1;
    }
}

void benchmark_lazywl_nf(string& s_rev, StreeOhleb<>& st, const size_type ntrials, size_type trial_length){
    size_type nt = 0;
    size_type i = 0;
    size_type k = s_rev.size() - 1;
    node_type v = st.root();

    while(nt++ < ntrials){
        for(i = 0; i < trial_length; i++)
            v = st.lazy_wl(v, s_rev[k--]);

        // SKIP finish completing the new node

        if(k < ntrials)
            k = s_rev.size() - 1;
    }
}

void benchmark_wl(string& s_rev, StreeOhleb<>& st, const size_type ntrials, size_type trial_length){
    size_type nt = 0;
    size_type i = 0;
    size_type k = s_rev.size() - 1;
    node_type v = st.root();

    while(nt++ < ntrials){
        for(i = 0; i < trial_length; i++)
            v = st.wl(v, s_rev[k--]);

        if(k < ntrials)
            k = s_rev.size() - 1;
    }
}

void benchmark_bstep(string& s_rev, StreeOhleb<>& st, const size_type ntrials, size_type trial_length){
    size_type nt = 0;
    size_type i = 0;
    size_type k = s_rev.size() - 1;
    node_type v = st.root();
    IInterval I = make_pair(v.i, v.j);

    while(nt++ < ntrials){
        for(i = 0; i < trial_length; i++){
            char c = s_rev[k--];
            int cc = st.csa.char2comp[c];
            I.first = st.csa.C[cc] + st.csa.bwt.rank(I.first, c);
            I.second = st.csa.C[cc] + st.csa.bwt.rank(I.second + 1, c) - 1;
        }
        if(k < ntrials)
            k = s_rev.size() - 1;
    }
}


int main(int argc, char **argv) {
    InputParser input(argc, argv);
    InputFlags flags(input);
    InputSpec S_fwd(input.getCmdOption("-s_path"));
    string s = S_fwd.load_s();

    Counter time_usage;
    StreeOhleb<> st;
    auto runs_start = timer::now();
    if(flags.load_stree){
        cerr << "loading the CST T(s) from " << S_fwd.s_fname + ".fwd.stree" << "... ";
        sdsl::load_from_file(st, S_fwd.s_fname + ".fwd.stree");
    } else {
        cerr << "building the CST T(s) of lentgth " << s.size() << "... ";
        sdsl::construct_im(st, s, 1);
    }
    auto runs_stop = timer::now();
    time_usage["dstruct"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    cerr << "DONE (" << time_usage["dstruct"] / 1000 << " seconds)" << endl;


    size_type ntrials = 10000;
    cout << "s_path,s,ntrials,nwlcalls,item,time_ms" << endl;
    for(size_type i = 1; i < 5; i = i + 1){
        cerr << ntrials << " of length : " << i << endl;
        auto start_time = timer::now();
        benchmark_wl(s, st, ntrials, i);
        time_usage["nonlazy"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
        cerr << " * nonlazy calls of length " << i << " took " << time_usage["nonlazy"] << " ms" << endl;

        start_time = timer::now();
        benchmark_lazywl(s, st, ntrials, i);
        time_usage["lazy"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
        cerr << " * lazy calls of length " << i << " took " << time_usage["lazy"] << " ms" << endl;

        start_time = timer::now();
        benchmark_lazywl_nf(s, st, ntrials, i);
        time_usage["lazy_nf"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
        cerr << " * lazy_nf calls of length " << i << " took " << time_usage["lazy_nf"] << " ms" << endl;

        start_time = timer::now();
        benchmark_bstep(s, st, ntrials, i);
        time_usage["bstep"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
        cerr << " * backward steps of length " << i << " took " << time_usage["bstep"] << " ms" << endl;

        for(auto item : time_usage){
            cout << S_fwd.s_fname << "," << s.size() << "," << ntrials << "," << i << "," << item.first << "," << item.second << endl;
        }
    }
    return 0;
}
